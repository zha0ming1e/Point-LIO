#include "li_initialization.h"
bool data_accum_finished = false, data_accum_start = false, online_calib_finish = false, refine_print = false;
int frame_num_init = 0;
double time_lag_IMU_wtr_lidar = 0.0, move_start_time = 0.0, online_calib_starts_time = 0.0; //, mean_acc_norm = 9.81;
double imu_first_time = 0.0;
bool lose_lid = false;
double timediff_imu_wrt_lidar = 0.0;
bool timediff_set_flg = false;
V3D gravity_lio = V3D::Zero();
mutex mtx_buffer;
sensor_msgs::Imu imu_last, imu_next;
// sensor_msgs::Imu::ConstPtr imu_last_ptr;
PointCloudXYZI::Ptr  ptr_con(new PointCloudXYZI());
double T1[MAXN], s_plot[MAXN], s_plot2[MAXN], s_plot3[MAXN], s_plot11[MAXN];

condition_variable sig_buffer;
int scan_count = 0;
int frame_ct = 0, wait_num = 0;
std::mutex m_time;
bool lidar_pushed = false, imu_pushed = false;
std::deque<PointCloudXYZI::Ptr>  lidar_buffer;
std::deque<double>               time_buffer;
std::deque<sensor_msgs::Imu::Ptr> imu_deque;
std::queue<std::vector<ObsPtr>> gnss_meas_buf;

void gnss_ephem_callback(const GnssEphemMsgConstPtr &ephem_msg)
{
    EphemPtr ephem = msg2ephem(ephem_msg);
    p_gnss->inputEphem(ephem);
}

void gnss_glo_ephem_callback(const GnssGloEphemMsgConstPtr &glo_ephem_msg)
{
    GloEphemPtr glo_ephem = msg2glo_ephem(glo_ephem_msg);
    p_gnss->inputEphem(glo_ephem);
}

void gnss_iono_params_callback(const StampedFloat64ArrayConstPtr &iono_msg)
{
    double ts = iono_msg->header.stamp.toSec();
    std::vector<double> iono_params;
    std::copy(iono_msg->data.begin(), iono_msg->data.end(), std::back_inserter(iono_params));
    assert(iono_params.size() == 8);
    p_gnss->inputIonoParams(ts, iono_params);
}

void rtk_pvt_callback(const GnssPVTSolnMsgConstPtr &groundt_pvt)
{
    double ts = time2sec(gst2time(groundt_pvt->time.week, groundt_pvt->time.tow));
    p_gnss->inputpvt(ts, groundt_pvt->latitude, groundt_pvt->longitude, groundt_pvt->altitude);
}

void rtk_lla_callback(const sensor_msgs::NavSatFixConstPtr &lla_msg)
{
    double ts = lla_msg->header.stamp.toSec();
    p_gnss->inputlla(ts, lla_msg->latitude, lla_msg->longitude, lla_msg->altitude);
}

void gnss_meas_callback(const GnssMeasMsgConstPtr &meas_msg)
{
    std::vector<ObsPtr> gnss_meas = msg2meas(meas_msg);
    
    latest_gnss_time = time2sec(gnss_meas[0]->time);
    // printf("gnss time: %f\n", latest_gnss_time);

    // cerr << "gnss ts is " << std::setprecision(20) << time2sec(gnss_meas[0]->time) << endl;
    if (!time_diff_valid)   return;

    mtx_buffer.lock();
    gnss_meas_buf.push(std::move(gnss_meas)); // ?
    mtx_buffer.unlock();
    sig_buffer.notify_all(); // notify_one()?
}

void local_trigger_info_callback(const ligo::LocalSensorExternalTriggerConstPtr &trigger_msg) // pps time sync
{
    std::lock_guard<std::mutex> lg(m_time);

    if (next_pulse_time_valid)
    {
        time_diff_gnss_local = next_pulse_time - trigger_msg->header.stamp.toSec();
        p_gnss->inputGNSSTimeDiff(time_diff_gnss_local);
        if (!time_diff_valid)       // just get calibrated
            std::cout << "time difference between GNSS and LI-Sensor got calibrated: "
                << std::setprecision(15) << time_diff_gnss_local << " s\n";
        time_diff_valid = true;
    }
}

void gnss_tp_info_callback(const GnssTimePulseInfoMsgConstPtr &tp_msg) // time stamp of GNSS signal
{
    gtime_t tp_time = gpst2time(tp_msg->time.week, tp_msg->time.tow);
    if (tp_msg->utc_based || tp_msg->time_sys == SYS_GLO)
        tp_time = utc2gpst(tp_time);
    else if (tp_msg->time_sys == SYS_GAL)
        tp_time = gst2time(tp_msg->time.week, tp_msg->time.tow);
    else if (tp_msg->time_sys == SYS_BDS)
        tp_time = bdt2time(tp_msg->time.week, tp_msg->time.tow);
    else if (tp_msg->time_sys == SYS_NONE)
    {
        std::cerr << "Unknown time system in GNSSTimePulseInfoMsg.\n";
        return;
    }
    double gnss_ts = time2sec(tp_time);

    std::lock_guard<std::mutex> lg(m_time);
    next_pulse_time = gnss_ts;
    next_pulse_time_valid = true;
}

void standard_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg) 
{
    mtx_buffer.lock();
    scan_count ++;
    double preprocess_start_time = omp_get_wtime();
    if (msg->header.stamp.toSec() < last_timestamp_lidar)
    {
        ROS_ERROR("lidar loop back, clear buffer");
        // lidar_buffer.shrink_to_fit();

        mtx_buffer.unlock();
        sig_buffer.notify_all();
        return;
    }

    last_timestamp_lidar = msg->header.stamp.toSec();
    // printf("check lidar time %f\n", last_timestamp_lidar);
    if (abs(last_timestamp_imu - last_timestamp_lidar) > 1.0 && !timediff_set_flg && !imu_deque.empty()) {
        timediff_set_flg = true;
        timediff_imu_wrt_lidar = last_timestamp_imu - last_timestamp_lidar;
        printf("Self sync IMU and LiDAR, HARD time lag is %.10lf \n \n", timediff_imu_wrt_lidar);
    }

    if ((lidar_type == VELO16 || lidar_type == OUST64 || lidar_type == HESAIxt32) && cut_frame_init) {
        deque<PointCloudXYZI::Ptr> ptr;
        deque<double> timestamp_lidar;
        p_pre->process_cut_frame_pcl2(msg, ptr, timestamp_lidar, cut_frame_num, scan_count);
        while (!ptr.empty() && !timestamp_lidar.empty()) {
            lidar_buffer.push_back(ptr.front());
            ptr.pop_front();
            time_buffer.push_back(timestamp_lidar.front() / double(1000));//unit:s
            timestamp_lidar.pop_front();
        }
    }
    else
    {
    PointCloudXYZI::Ptr  ptr(new PointCloudXYZI(20000,1));
    p_pre->process(msg, ptr);
    if (con_frame)
    {
        if (frame_ct == 0)
        {
            time_con = last_timestamp_lidar; //msg->header.stamp.toSec();
        }
        if (frame_ct < 10)
        {
            for (int i = 0; i < ptr->size(); i++)
            {
                ptr->points[i].curvature += (last_timestamp_lidar - time_con) * 1000;
                ptr_con->push_back(ptr->points[i]);
            }
            frame_ct ++;
        }
        else
        {
            PointCloudXYZI::Ptr  ptr_con_i(new PointCloudXYZI(10000,1));
            // cout << "ptr div num:" << ptr_div->size() << endl;
            *ptr_con_i = *ptr_con;
            lidar_buffer.push_back(ptr_con_i);
            double time_con_i = time_con;
            time_buffer.push_back(time_con_i);
            ptr_con->clear();
            frame_ct = 0;
        }
    }
    else
    { 
        lidar_buffer.emplace_back(ptr);
        time_buffer.emplace_back(msg->header.stamp.toSec());
    }
    }
    s_plot11[scan_count] = omp_get_wtime() - preprocess_start_time;
    mtx_buffer.unlock();
    sig_buffer.notify_all();
}

void livox_pcl_cbk(const livox_ros_driver::CustomMsg::ConstPtr &msg) 
{
    mtx_buffer.lock();
    double preprocess_start_time = omp_get_wtime();
    scan_count ++;
    if (msg->header.stamp.toSec() < last_timestamp_lidar)
    {
        ROS_ERROR("lidar loop back, clear buffer");

        mtx_buffer.unlock();
        sig_buffer.notify_all();
        return;
        // lidar_buffer.shrink_to_fit();
    }

    last_timestamp_lidar = msg->header.stamp.toSec();    
    if (abs(last_timestamp_imu - last_timestamp_lidar) > 1.0 && !timediff_set_flg && !imu_deque.empty()) {
        timediff_set_flg = true;
        timediff_imu_wrt_lidar = last_timestamp_imu - last_timestamp_lidar;
        printf("Self sync IMU and LiDAR, HARD time lag is %.10lf \n \n", timediff_imu_wrt_lidar);
    }

    if (cut_frame_init) {
        deque<PointCloudXYZI::Ptr> ptr;
        deque<double> timestamp_lidar;
        p_pre->process_cut_frame_livox(msg, ptr, timestamp_lidar, cut_frame_num, scan_count);

        while (!ptr.empty() && !timestamp_lidar.empty()) {
            lidar_buffer.push_back(ptr.front());
            ptr.pop_front();
            time_buffer.push_back(timestamp_lidar.front() / double(1000));//unit:s
            timestamp_lidar.pop_front();
        }
    }
    else
    {
    PointCloudXYZI::Ptr  ptr(new PointCloudXYZI(10000,1));
    p_pre->process(msg, ptr); 
    if (con_frame)
    {
        if (frame_ct == 0)
        {
            time_con = last_timestamp_lidar; //msg->header.stamp.toSec();
        }
        if (frame_ct < 10)
        {
            for (int i = 0; i < ptr->size(); i++)
            {
                ptr->points[i].curvature += (last_timestamp_lidar - time_con) * 1000;
                ptr_con->push_back(ptr->points[i]);
            }
            frame_ct ++;
        }
        else
        {
            PointCloudXYZI::Ptr  ptr_con_i(new PointCloudXYZI(10000,1));
            // cout << "ptr div num:" << ptr_div->size() << endl;
            *ptr_con_i = *ptr_con;
            double time_con_i = time_con;
            lidar_buffer.push_back(ptr_con_i);
            time_buffer.push_back(time_con_i);
            ptr_con->clear();
            frame_ct = 0;
        }
    }
    else
    {
        lidar_buffer.emplace_back(ptr);
        time_buffer.emplace_back(msg->header.stamp.toSec());
    }
    }
    s_plot11[scan_count] = omp_get_wtime() - preprocess_start_time;
    mtx_buffer.unlock();
    sig_buffer.notify_all();
}

void imu_cbk(const sensor_msgs::Imu::ConstPtr &msg_in) 
{
    mtx_buffer.lock();

    // publish_count ++;
    sensor_msgs::Imu::Ptr msg(new sensor_msgs::Imu(*msg_in));

    msg->header.stamp = ros::Time().fromSec(msg->header.stamp.toSec() - timediff_imu_wrt_lidar - time_lag_IMU_wtr_lidar);

    double timestamp = msg->header.stamp.toSec();
    // printf("time_diff%f, %f, %f\n", last_timestamp_imu - timestamp, last_timestamp_imu, timestamp);

    if (timestamp < last_timestamp_imu)
    {
        ROS_ERROR("imu loop back, clear deque");
        // imu_deque.shrink_to_fit();
        // cout << "check time:" << timestamp << ";" << last_timestamp_imu << endl;
        // printf("time_diff%f, %f, %f\n", last_timestamp_imu - timestamp, last_timestamp_imu, timestamp);
        if (p_imu->UseLIInit)
        {
            Init_LI->IMU_buffer_clear();
        }
        mtx_buffer.unlock();
        sig_buffer.notify_all();
        return;
    }
    // push all IMU meas into Init_LI
    if (!p_imu->LI_init_done && !data_accum_finished)
        Init_LI->push_ALL_IMU_CalibState(msg, acc_norm); // mean_acc_norm);

    imu_deque.emplace_back(msg);
    last_timestamp_imu = timestamp;
    mtx_buffer.unlock();
    sig_buffer.notify_all();
}

bool sync_packages(MeasureGroup &meas, queue<std::vector<ObsPtr>> &gnss_msg)
{
    if (nolidar)
    {
        if (imu_deque.empty())
        {
            return false;
        }
        
        imu_first_time = imu_deque.front()->header.stamp.toSec(); // 

        if (latest_gnss_time < gnss_local_time_diff + imu_first_time + lidar_time_inte)
        {
            return false;
        }

        if (last_timestamp_imu < imu_first_time + lidar_time_inte)
        {
            return false;
        }

        if (is_first_gnss)
        {
            if (gnss_meas_buf.empty())
            {
                // imu_pushed = false;
                return false;
            }
            else
            {
                double imu_time = imu_deque.front()->header.stamp.toSec();
                double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
                if (last_timestamp_imu < front_gnss_ts)
                {
                    return false;
                }
                while (front_gnss_ts - imu_time > gnss_local_time_diff) // wrong
                {
                    imu_last = *(imu_deque.front());
                    imu_deque.pop_front();
                    if(imu_deque.empty()) break;
                    imu_time = imu_deque.front()->header.stamp.toSec(); // can be changed
                    imu_next = *(imu_deque.front());
                }
                // else
                {
                    is_first_gnss = false;
                }
            }
        }

        if (imu_deque.empty())
        {
            cout << "could not be here" << endl;
            return false;
        }

        if (!imu_pushed)
        { 
            double imu_time = imu_deque.front()->header.stamp.toSec();
            // imu_first_time = imu_time;

            double imu_last_time = imu_deque.back()->header.stamp.toSec();
            if (imu_last_time - imu_first_time < lidar_time_inte)
            {
                return false;
            }
            /*** push imu data, and pop from imu buffer ***/
            if (p_imu->imu_need_init_)
            {
                imu_next = *(imu_deque.front());
                meas.imu.shrink_to_fit();
                while (imu_time - imu_first_time < lidar_time_inte)
                {
                    meas.imu.emplace_back(imu_deque.front());
                    imu_last = imu_next;
                    imu_deque.pop_front();
                    if(imu_deque.empty()) break;
                    imu_time = imu_deque.front()->header.stamp.toSec(); // can be changed
                    imu_next = *(imu_deque.front());
                }
                if (!gnss_meas_buf.empty())
                {
                    double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
                    while (front_gnss_ts < imu_first_time + lidar_time_inte + gnss_local_time_diff)
                    {
                        gnss_meas_buf.pop();
                        if(gnss_meas_buf.empty()) break;
                        front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
                    }
                    if (meas.imu.empty())
                    {
                        return false;
                    }
                }
            }
            imu_pushed = true;
        }

        if (GNSS_ENABLE)
        {
            if (!gnss_meas_buf.empty()) // or can wait for a short time?
            {
                // double back_gnss_ts = time2sec(gnss_meas_buf.back()[0]->time);
                
                // if (back_gnss_ts - imu_first_time < gnss_local_time_diff + lidar_time_inte)
                // {
                //     return false;
                // }
                double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
                while (front_gnss_ts - imu_first_time < gnss_local_time_diff + lidar_time_inte) 
                {
                    gnss_msg.push(gnss_meas_buf.front());
                    gnss_meas_buf.pop();
                    if (gnss_meas_buf.empty()) break;
                    front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
                }

                if (!gnss_msg.empty())
                {
                    imu_pushed = false;
                    return true;
                }
            }

            if (gnss_meas_buf.empty())
            {
                wait_num ++;
                if (wait_num > 5) 
                {
                    wait_num = 0;
                }
                else
                {
                    return false;
                }
            }
        }
        imu_pushed = false;
        return true;
    }
    else
    {
    if (!imu_en)
    {
        if (!lidar_buffer.empty())
        {
            if (!lidar_pushed)
            {
                meas.lidar = lidar_buffer.front();
                meas.lidar_beg_time = time_buffer.front();
                lose_lid = false;
                if(meas.lidar->points.size() < 1) 
                {
                    cout << "lose lidar" << std::endl;
                    // return false;
                    lose_lid = true;
                }
                else
                {
                    double end_time = meas.lidar->points.back().curvature;
                    for (auto pt: meas.lidar->points)
                    {
                        if (pt.curvature > end_time)
                        {
                            end_time = pt.curvature;
                        }
                    }
                    lidar_end_time = meas.lidar_beg_time + end_time / double(1000);
                    meas.lidar_last_time = lidar_end_time;
                }
                lidar_pushed = true;
            }
            
            if (GNSS_ENABLE)
            {
                if (!gnss_meas_buf.empty()) // or can wait for a short time?
                {
                    // double back_gnss_ts = time2sec(gnss_meas_buf.back()[0]->time);
                    // if (!lose_lid && back_gnss_ts < lidar_end_time)
                    // {
                    //     return false;
                    // }
                    // if (lose_lid && back_gnss_ts < meas.lidar_beg_time + lidar_time_inte)
                    // {
                    //     return false;
                    // }
                    double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
                    while (front_gnss_ts < meas.lidar_beg_time + gnss_local_time_diff) // 0.05
                    {
                        ROS_WARN("throw gnss, only should happen at the beginning 542");
                        gnss_meas_buf.pop();
                        if (gnss_meas_buf.empty()) break;
                        front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
                    }
                    while (((!lose_lid && front_gnss_ts < lidar_end_time + gnss_local_time_diff) || (lose_lid && front_gnss_ts < meas.lidar_beg_time + gnss_local_time_diff + lidar_time_inte) ))
                    {
                        gnss_msg.push(gnss_meas_buf.front());
                        gnss_meas_buf.pop();
                        if (gnss_meas_buf.empty()) break;
                        front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
                    }
                    if (!gnss_msg.empty())
                    {
                        time_buffer.pop_front();
                        lidar_buffer.pop_front();
                        lidar_pushed = false;
                        return true;
                    }
                }

                if (gnss_meas_buf.empty())
                {
                    wait_num ++;
                    if (wait_num > 5) 
                    {
                        wait_num = 0;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
            time_buffer.pop_front();
            lidar_buffer.pop_front();
            lidar_pushed = false;
            if (!lose_lid)
            {
                return true;
            }
            else
            {
                return false;
            }
        }        
        return false;
    }

    if (lidar_buffer.empty() || imu_deque.empty())
    {
        return false;
    }
    /*** push a lidar scan ***/
    if(!lidar_pushed)
    {
        lose_lid = false;
        meas.lidar = lidar_buffer.front();
        meas.lidar_beg_time = time_buffer.front();
        if(meas.lidar->points.size() < 1) 
        {
            cout << "lose lidar" << endl;
            lose_lid = true;
            lidar_pushed = true;
            // lidar_buffer.pop_front();
            // time_buffer.pop_front();
            // return false;
        }
        else
        {
            double end_time = meas.lidar->points.back().curvature;
            for (auto pt: meas.lidar->points)
            {
                if (pt.curvature > end_time)
                {
                    end_time = pt.curvature;
                }
            }
            lidar_end_time = meas.lidar_beg_time + end_time / double(1000);
            // cout << "check time lidar:" << end_time << endl;
            meas.lidar_last_time = lidar_end_time;
            lidar_pushed = true;
        }
    }

    if (!lose_lid && last_timestamp_imu < lidar_end_time)
    {
        return false;
    }
    if (lose_lid && last_timestamp_imu < meas.lidar_beg_time + lidar_time_inte)
    {
        return false;
    }

    if (is_first_gnss)
    {
        if (GNSS_ENABLE)
        {
            if (gnss_meas_buf.empty() && lose_lid)
            {
                lidar_pushed = false;

                lidar_buffer.pop_front();
                time_buffer.pop_front();
                return false;
            }
            else
            {
                if (!lose_lid)
                {
                    double imu_time = imu_deque.front()->header.stamp.toSec();
                    imu_next = *(imu_deque.front());
                    // meas.imu.shrink_to_fit();
                    while (imu_time < meas.lidar_beg_time)
                    {
                        imu_last = imu_next;
                        imu_deque.pop_front();
                        if(imu_deque.empty()) break;
                        imu_time = imu_deque.front()->header.stamp.toSec(); // can be changed
                        imu_next = *(imu_deque.front());
                    }
                }
                else if (!gnss_meas_buf.empty())
                {
                    double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time

                    double imu_time = imu_deque.front()->header.stamp.toSec();
                    imu_next = *(imu_deque.front());
                    // meas.imu.shrink_to_fit();
                    if (front_gnss_ts > meas.lidar_beg_time + lidar_time_inte + gnss_local_time_diff)
                    {
                        while (imu_time < meas.lidar_beg_time + lidar_time_inte)
                        {
                            imu_last = imu_next;
                            imu_deque.pop_front();
                            if(imu_deque.empty()) break;
                            imu_time = imu_deque.front()->header.stamp.toSec(); // can be changed
                            imu_next = *(imu_deque.front());
                        }
                        lidar_pushed = false;
                        lidar_buffer.pop_front();
                        time_buffer.pop_front();
                        return false;
                    }
                    else
                    {
                        while (imu_time < front_gnss_ts - gnss_local_time_diff)
                        {
                            imu_deque.pop_front();
                            if(imu_deque.empty()) break;
                            imu_time = imu_deque.front()->header.stamp.toSec(); // can be changed
                            imu_next = *(imu_deque.front());
                        }
                    }
                }
                {
                    if (lose_lid)
                    {
                        if (imu_deque.empty())
                        {
                            lidar_pushed = false;
                            lidar_buffer.pop_front();
                            time_buffer.pop_front();
                            return false;
                        }
                        else if (imu_deque.front()->header.stamp.toSec() > meas.lidar_beg_time + lidar_time_inte)
                        {
                            lidar_pushed = false;
                            lidar_buffer.pop_front();
                            time_buffer.pop_front();
                            is_first_gnss = false;
                            return false;
                        }
                    }
                    is_first_gnss = false;
                }
            }
        }
        else
        {
            if (lose_lid)
            {
                lidar_pushed = false;

                lidar_buffer.pop_front();
                time_buffer.pop_front();
                return false;
            }
            else
            {
                // if (!lose_lid)
                {
                    double imu_time = imu_deque.front()->header.stamp.toSec();
                    imu_next = *(imu_deque.front());
                    // meas.imu.shrink_to_fit();
                    while (imu_time < meas.lidar_beg_time)
                    {
                        imu_last = imu_next;
                        imu_deque.pop_front();
                        if(imu_deque.empty()) break;
                        imu_time = imu_deque.front()->header.stamp.toSec(); // can be changed
                        imu_next = *(imu_deque.front());
                    }
                }
                {
                    is_first_gnss = false;
                }
            }
        }
    }

    // if (!lose_lid && last_timestamp_imu < lidar_end_time)
    // {
    //     return false;
    // }
    // if (lose_lid && last_timestamp_imu < meas.lidar_beg_time + lidar_time_inte)
    // {
    //     return false;
    // }

    if (!lose_lid && !imu_pushed)
    { 
        /*** push imu data, and pop from imu buffer ***/
        if (p_imu->imu_need_init_)
        {
            double imu_time = imu_deque.front()->header.stamp.toSec();
            imu_next = *(imu_deque.front());
            meas.imu.shrink_to_fit();
            while (imu_time < lidar_end_time)
            {
                meas.imu.emplace_back(imu_deque.front());
                imu_last = imu_next;
                imu_deque.pop_front();
                if(imu_deque.empty()) break;
                imu_time = imu_deque.front()->header.stamp.toSec(); // can be changed
                imu_next = *(imu_deque.front());
            }
            if (GNSS_ENABLE)
            {
                if (!gnss_meas_buf.empty())
                {
                    double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
                    while (front_gnss_ts < lidar_end_time + gnss_local_time_diff)
                    {
                        gnss_meas_buf.pop();
                        if(gnss_meas_buf.empty()) break;
                        front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
                    }
                }
            }
        }
        imu_pushed = true;
    }

    if (lose_lid && !imu_pushed)
    { 
        /*** push imu data, and pop from imu buffer ***/
        if (p_imu->imu_need_init_)
        {
            double imu_time = imu_deque.front()->header.stamp.toSec();
            meas.imu.shrink_to_fit();

            imu_next = *(imu_deque.front());
            while (imu_time < meas.lidar_beg_time + lidar_time_inte)
            {
                meas.imu.emplace_back(imu_deque.front());
                imu_last = imu_next;
                imu_deque.pop_front();
                if(imu_deque.empty()) break;
                imu_time = imu_deque.front()->header.stamp.toSec(); // can be changed
                imu_next = *(imu_deque.front());
            }

            if (GNSS_ENABLE)
            {
                if (!gnss_meas_buf.empty())
                {
                    double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
                    while (front_gnss_ts < meas.lidar_beg_time + lidar_time_inte + gnss_local_time_diff)
                    {
                        gnss_meas_buf.pop();
                        if(gnss_meas_buf.empty()) break;
                        front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
                    }
                }
            }
        }
        imu_pushed = true;
    }

    if (GNSS_ENABLE)
    {
        if (!gnss_meas_buf.empty()) // or can wait for a short time?
        {
            // double back_gnss_ts = time2sec(gnss_meas_buf.back()[0]->time);
            // if (!lose_lid && back_gnss_ts < lidar_end_time)
            // {
            //     return false;
            // }
            // if (lose_lid && back_gnss_ts < meas.lidar_beg_time + lidar_time_inte)
            // {
            //     return false;
            // }
            double front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time); // take time
            while (((!lose_lid && front_gnss_ts < lidar_end_time + gnss_local_time_diff) || (lose_lid && front_gnss_ts < meas.lidar_beg_time + gnss_local_time_diff + lidar_time_inte) ))
            {
                gnss_msg.push(gnss_meas_buf.front());
                gnss_meas_buf.pop();
                if (gnss_meas_buf.empty()) break;
                front_gnss_ts = time2sec(gnss_meas_buf.front()[0]->time);
            }
            if (!gnss_msg.empty())
            {
                time_buffer.pop_front();
                lidar_buffer.pop_front();
                lidar_pushed = false;
                imu_pushed = false;
                return true;
            }
        }

        if (gnss_meas_buf.empty())
        {
            wait_num ++;
            if (wait_num > 5) 
            {
                wait_num = 0;
            }
            else
            {
                return false;
            }
        }
    }

    lidar_buffer.pop_front();
    time_buffer.pop_front();
    lidar_pushed = false;
    imu_pushed = false;
    return true;
    // if (!lose_lid)
    // {
    //     return true;
    // }
    // else
    // {
    //     return false;
    // }
    }
}

void LI_Init_set()
{
    //Transfer to PointLIO
    // state.offset_R_L_I = Init_LI->get_R_LI();
    // state.offset_T_L_I = Init_LI->get_T_LI();
    Lidar_R_wrt_IMU = Init_LI->get_R_LI();
    Lidar_T_wrt_IMU = Init_LI->get_T_LI();
    if (use_imu_as_input)
    {
        kf_input.x_.pos = -(p_imu->state_LI_Init.rot.normalized().toRotationMatrix() * Lidar_R_wrt_IMU.transpose() * Lidar_T_wrt_IMU +
                        p_imu->state_LI_Init.pos); //Body frame is IMU frame in Point-LIO mode
        kf_input.x_.rot = p_imu->state_LI_Init.rot.normalized().toRotationMatrix() * Lidar_R_wrt_IMU.transpose();
        gravity_lio = Init_LI->get_Grav_L0();
        kf_input.x_.gravity = Init_LI->get_Grav_L0();
        kf_input.x_.bg = Init_LI->get_gyro_bias();
        kf_input.x_.ba = Init_LI->get_acc_bias();
    }
    else
    {
        kf_output.x_.pos = -1 * (p_imu->state_LI_Init.rot.normalized().toRotationMatrix() * Lidar_R_wrt_IMU.transpose() * Lidar_T_wrt_IMU +
                        p_imu->state_LI_Init.pos); //Body frame is IMU frame in Point-LIO mode
        kf_output.x_.rot = p_imu->state_LI_Init.rot.normalized().toRotationMatrix() * Lidar_R_wrt_IMU.transpose();
        gravity_lio = Init_LI->get_Grav_L0();
        kf_output.x_.gravity = Init_LI->get_Grav_L0();
        kf_output.x_.bg = Init_LI->get_gyro_bias();
        kf_output.x_.ba = Init_LI->get_acc_bias();
        kf_output.x_.acc = - kf_output.x_.rot.normalized().toRotationMatrix().transpose() * gravity_lio;
    }
    {
        // init_map = true;
        // delete ikdtree.Root_Node;
        // ikdtree.Root_Node = nullptr;
    }
    if (GNSS_ENABLE)
    {
        p_imu->Set_init(Rot_gnss_init); // gravity_lio, 
        p_gnss->Rot_gnss_init = Rot_gnss_init;
        cout << "check Rot init:" << Rot_gnss_init << endl;
        // delete ikdtree.Root_Node;
        // ikdtree.Root_Node = nullptr;

    }
    if (lidar_type != AVIA)
        cut_frame_num = 2;

    time_lag_IMU_wtr_lidar = Init_LI->get_total_time_lag(); //Compensate IMU's time in the buffer
    for (int i = 0; i < imu_deque.size(); i++) {
        imu_deque[i]->header.stamp = ros::Time().fromSec(imu_deque[i]->header.stamp.toSec()- time_lag_IMU_wtr_lidar);
    }
    // last_timestamp_imu2 = last_timestamp_imu2 - time_lag_IMU_wtr_lidar; // He

    p_imu->LI_init_done = true;
    p_imu->imu_need_init_ = false;
    cut_frame_init = false;
}

