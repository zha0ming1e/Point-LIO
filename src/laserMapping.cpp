// #include <so3_math.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <visualization_msgs/Marker.h>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/io/pcd_io.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>
#include "li_initialization.h"
// #include "Estimator.h"
#include <malloc.h>
// #include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include "chi-square.h"
// #include "matplotlibcpp.h"
// #include "LI_init/LI_init.h"
// #include <ros/console.h>


#define PUBFRAME_PERIOD     (20)

const float MOV_THRESHOLD = 1.5f;

string root_dir = ROOT_DIR;

int time_log_counter = 0, publish_count = 0;

bool init_map = false, flg_first_scan = true;
std::vector<ObsPtr> gnss_cur;

// Time Log Variables
double match_time = 0, solve_time = 0, propag_time = 0, update_time = 0;

bool  flg_reset = false, flg_exit = false;

vector<BoxPointType> cub_needrm;

//surf feature in map
PointCloudXYZI::Ptr feats_undistort(new PointCloudXYZI());
PointCloudXYZI::Ptr feats_down_body_space(new PointCloudXYZI());
PointCloudXYZI::Ptr init_feats_world(new PointCloudXYZI());

pcl::VoxelGrid<PointType> downSizeFilterSurf;
pcl::VoxelGrid<PointType> downSizeFilterMap;

V3D euler_cur;

nav_msgs::Path path;
nav_msgs::Odometry odomAftMapped;
geometry_msgs::PoseStamped msg_body_pose;

void SigHandle(int sig)
{
    flg_exit = true;
    ROS_WARN("catch sig %d", sig);
    sig_buffer.notify_all();
}

inline void dump_lio_state_to_log(FILE *fp)  
{
    V3D rot_ang;
    if (!use_imu_as_input)
    {
        rot_ang = SO3ToEuler(kf_output.x_.rot);
    }
    else
    {
        rot_ang = SO3ToEuler(kf_input.x_.rot);
    }
    
    fprintf(fp, "%lf ", Measures.lidar_beg_time - first_lidar_time);
    fprintf(fp, "%lf %lf %lf ", rot_ang(0), rot_ang(1), rot_ang(2));                   // Angle
    if (use_imu_as_input)
    {
        fprintf(fp, "%lf %lf %lf ", kf_input.x_.pos(0), kf_input.x_.pos(1), kf_input.x_.pos(2)); // Pos  
        fprintf(fp, "%lf %lf %lf ", 0.0, 0.0, 0.0);                                        // omega  
        fprintf(fp, "%lf %lf %lf ", kf_input.x_.vel(0), kf_input.x_.vel(1), kf_input.x_.vel(2)); // Vel  
        fprintf(fp, "%lf %lf %lf ", 0.0, 0.0, 0.0);                                        // Acc  
        fprintf(fp, "%lf %lf %lf ", kf_input.x_.bg(0), kf_input.x_.bg(1), kf_input.x_.bg(2));    // Bias_g  
        fprintf(fp, "%lf %lf %lf ", kf_input.x_.ba(0), kf_input.x_.ba(1), kf_input.x_.ba(2));    // Bias_a  
        fprintf(fp, "%lf %lf %lf ", kf_input.x_.gravity(0), kf_input.x_.gravity(1), kf_input.x_.gravity(2)); // Bias_a  
    }
    else
    {
        fprintf(fp, "%lf %lf %lf ", kf_output.x_.pos(0), kf_output.x_.pos(1), kf_output.x_.pos(2)); // Pos  
        fprintf(fp, "%lf %lf %lf ", 0.0, 0.0, 0.0);                                        // omega  
        fprintf(fp, "%lf %lf %lf ", kf_output.x_.vel(0), kf_output.x_.vel(1), kf_output.x_.vel(2)); // Vel  
        fprintf(fp, "%lf %lf %lf ", 0.0, 0.0, 0.0);                                        // Acc  
        fprintf(fp, "%lf %lf %lf ", kf_output.x_.bg(0), kf_output.x_.bg(1), kf_output.x_.bg(2));    // Bias_g  
        fprintf(fp, "%lf %lf %lf ", kf_output.x_.ba(0), kf_output.x_.ba(1), kf_output.x_.ba(2));    // Bias_a  
        fprintf(fp, "%lf %lf %lf ", kf_output.x_.gravity(0), kf_output.x_.gravity(1), kf_output.x_.gravity(2)); // Bias_a  
    }
    fprintf(fp, "\r\n");  
    fflush(fp);
}

void pointBodyLidarToIMU(PointType const * const pi, PointType * const po)
{
    V3D p_body_lidar(pi->x, pi->y, pi->z);
    V3D p_body_imu;
    // if (extrinsic_est_en)
    // {
    //     if (!use_imu_as_input)
    //     {
    //         p_body_imu = kf_output.x_.offset_R_L_I * p_body_lidar + kf_output.x_.offset_T_L_I;
    //     }
    //     else
    //     {
    //         p_body_imu = kf_input.x_.offset_R_L_I * p_body_lidar + kf_input.x_.offset_T_L_I;
    //     }
    // }
    // else
    {
        p_body_imu = Lidar_R_wrt_IMU * p_body_lidar + Lidar_T_wrt_IMU;
    }
    po->x = p_body_imu(0);
    po->y = p_body_imu(1);
    po->z = p_body_imu(2);
    po->intensity = pi->intensity;
}

int points_cache_size = 0;

void points_cache_collect() // seems for debug
{
    PointVector points_history;
    ikdtree.acquire_removed_points(points_history);
    points_cache_size = points_history.size();
}

BoxPointType LocalMap_Points;
bool Localmap_Initialized = false;
void lasermap_fov_segment()
{
    cub_needrm.shrink_to_fit();

    V3D pos_LiD;
    if (use_imu_as_input)
    {
        state_in = kf_input.x_;
        pos_LiD = state_in.pos + state_in.rot.normalized() * Lidar_T_wrt_IMU;
    }
    else
    {
        state_out = kf_output.x_;
        pos_LiD = state_out.pos + state_out.rot.normalized() * Lidar_T_wrt_IMU;
    }
    if (!Localmap_Initialized){
        for (int i = 0; i < 3; i++){
            LocalMap_Points.vertex_min[i] = pos_LiD(i) - cube_len / 2.0;
            LocalMap_Points.vertex_max[i] = pos_LiD(i) + cube_len / 2.0;
        }
        Localmap_Initialized = true;
        return;
    }
    float dist_to_map_edge[3][2];
    bool need_move = false;
    for (int i = 0; i < 3; i++){
        dist_to_map_edge[i][0] = fabs(pos_LiD(i) - LocalMap_Points.vertex_min[i]);
        dist_to_map_edge[i][1] = fabs(pos_LiD(i) - LocalMap_Points.vertex_max[i]);
        if (dist_to_map_edge[i][0] <= MOV_THRESHOLD * DET_RANGE || dist_to_map_edge[i][1] <= MOV_THRESHOLD * DET_RANGE) need_move = true;
    }
    if (!need_move) return;
    BoxPointType New_LocalMap_Points, tmp_boxpoints;
    New_LocalMap_Points = LocalMap_Points;
    float mov_dist = max((cube_len - 2.0 * MOV_THRESHOLD * DET_RANGE) * 0.5 * 0.9, double(DET_RANGE * (MOV_THRESHOLD -1)));
    for (int i = 0; i < 3; i++){
        tmp_boxpoints = LocalMap_Points;
        if (dist_to_map_edge[i][0] <= MOV_THRESHOLD * DET_RANGE){
            New_LocalMap_Points.vertex_max[i] -= mov_dist;
            New_LocalMap_Points.vertex_min[i] -= mov_dist;
            tmp_boxpoints.vertex_min[i] = LocalMap_Points.vertex_max[i] - mov_dist;
            cub_needrm.emplace_back(tmp_boxpoints);
        } else if (dist_to_map_edge[i][1] <= MOV_THRESHOLD * DET_RANGE){
            New_LocalMap_Points.vertex_max[i] += mov_dist;
            New_LocalMap_Points.vertex_min[i] += mov_dist;
            tmp_boxpoints.vertex_max[i] = LocalMap_Points.vertex_min[i] + mov_dist;
            cub_needrm.emplace_back(tmp_boxpoints);
        }
    }
    LocalMap_Points = New_LocalMap_Points;

    points_cache_collect();
}

int process_increments = 0;
void map_incremental()
{
    PointVector PointToAdd;
    PointVector PointNoNeedDownsample;
    PointToAdd.reserve(feats_down_size);
    PointNoNeedDownsample.reserve(feats_down_size);
    
        for(int i = 0; i < feats_down_size; i++)
        {
            if (p_imu->UseLIInit && !p_imu->LI_init_done)
            {
                p_imu->pointBodyToWorld_li_init(&(feats_down_body->points[i]), &(feats_down_world->points[i]));
            }
            /* No points found within the given threshold of nearest search*/
            if (Nearest_Points[i].empty()){
                
                PointNoNeedDownsample.emplace_back(feats_down_world->points[i]);
                continue;          
            }      
            /* decide if need add to map */
            
            if (!Nearest_Points[i].empty())
            {
                const PointVector &points_near = Nearest_Points[i];
                bool need_add = true;
                BoxPointType Box_of_Point;
                PointType downsample_result, mid_point; 
                mid_point.x = floor(feats_down_world->points[i].x/filter_size_map_min)*filter_size_map_min + 0.5 * filter_size_map_min;
                mid_point.y = floor(feats_down_world->points[i].y/filter_size_map_min)*filter_size_map_min + 0.5 * filter_size_map_min;
                mid_point.z = floor(feats_down_world->points[i].z/filter_size_map_min)*filter_size_map_min + 0.5 * filter_size_map_min;
                /* If the nearest points is definitely outside the downsample box */
                if (fabs(points_near[0].x - mid_point.x) > 1.732 * filter_size_map_min || fabs(points_near[0].y - mid_point.y) > 1.732 * filter_size_map_min || fabs(points_near[0].z - mid_point.z) > 1.732 * filter_size_map_min){
                    PointNoNeedDownsample.emplace_back(feats_down_world->points[i]);
                    continue;
                }
                /* Check if there is a point already in the downsample box and closer to the center point */
                float dist  = calc_dist<float>(feats_down_world->points[i],mid_point);
                for (int readd_i = 0; readd_i < NUM_MATCH_POINTS; readd_i ++)
                {
                    if (points_near.size() < NUM_MATCH_POINTS) break;
                    /* Those points which are outside the downsample box should not be considered. */
                    if (fabs(points_near[readd_i].x - mid_point.x) > 0.5 * filter_size_map_min || fabs(points_near[readd_i].y - mid_point.y) > 0.5 * filter_size_map_min || fabs(points_near[readd_i].z - mid_point.z) > 0.5 * filter_size_map_min) {
                        continue;                    
                    }
                    if (calc_dist<float>(points_near[readd_i], mid_point) < dist)
                    {
                        need_add = false;
                        break;
                    }
                }
                if (need_add) PointToAdd.emplace_back(feats_down_world->points[i]);
            }
            else
            {
                PointToAdd.emplace_back(feats_down_world->points[i]);
            }
        }
    
    ikdtree.Add_Points(PointNoNeedDownsample, false);
}

void publish_init_kdtree(const ros::Publisher & pubLaserCloudFullRes)
{
    int size_init_ikdtree = ikdtree.size();
    PointCloudXYZI::Ptr   laserCloudInit(new PointCloudXYZI(size_init_ikdtree, 1));

    sensor_msgs::PointCloud2 laserCloudmsg;
    PointVector ().swap(ikdtree.PCL_Storage);
    ikdtree.flatten(ikdtree.Root_Node, ikdtree.PCL_Storage, NOT_RECORD);
                
    laserCloudInit->points = ikdtree.PCL_Storage;
    pcl::toROSMsg(*laserCloudInit, laserCloudmsg);
        
    laserCloudmsg.header.stamp = ros::Time().fromSec(lidar_end_time);
    laserCloudmsg.header.frame_id = "camera_init";
    pubLaserCloudFullRes.publish(laserCloudmsg);

}

PointCloudXYZI::Ptr pcl_wait_pub(new PointCloudXYZI(500000, 1));
PointCloudXYZI::Ptr pcl_wait_save(new PointCloudXYZI());
void publish_frame_world(const ros::Publisher & pubLaserCloudFullRes)
{
    if (scan_pub_en)
    {
        PointCloudXYZI::Ptr laserCloudFullRes(feats_down_body);
        int size = laserCloudFullRes->points.size() / 3;

        PointCloudXYZI::Ptr   laserCloudWorld(new PointCloudXYZI(size, 1));
        
        for (int i = 0; i < size; i++)
        {
            // if (i % 3 == 0)
            // {
            laserCloudWorld->points[i].x = feats_down_world->points[i*3].x;
            laserCloudWorld->points[i].y = feats_down_world->points[i*3].y;
            laserCloudWorld->points[i].z = feats_down_world->points[i*3].z;
            laserCloudWorld->points[i].intensity = feats_down_world->points[i*3].intensity; // feats_down_world->points[i].y; // 
            // }
        }
        sensor_msgs::PointCloud2 laserCloudmsg;
        pcl::toROSMsg(*laserCloudWorld, laserCloudmsg);
        
        laserCloudmsg.header.stamp = ros::Time().fromSec(lidar_end_time);
        laserCloudmsg.header.frame_id = "camera_init";
        pubLaserCloudFullRes.publish(laserCloudmsg);
        publish_count -= PUBFRAME_PERIOD;
    }
    
    /**************** save map ****************/
    /* 1. make sure you have enough memories
    /* 2. noted that pcd save will influence the real-time performences **/
    if (pcd_save_en)
    {
        int size = feats_down_world->points.size();
        PointCloudXYZI::Ptr   laserCloudWorld(new PointCloudXYZI(size, 1));

        for (int i = 0; i < size; i++)
        {
            laserCloudWorld->points[i].x = feats_down_world->points[i].x;
            laserCloudWorld->points[i].y = feats_down_world->points[i].y;
            laserCloudWorld->points[i].z = feats_down_world->points[i].z;
            laserCloudWorld->points[i].intensity = feats_down_world->points[i].intensity;
        }

        *pcl_wait_save += *laserCloudWorld;

        static int scan_wait_num = 0;
        scan_wait_num ++;
        if (pcl_wait_save->size() > 0 && pcd_save_interval > 0  && scan_wait_num >= pcd_save_interval)
        {
            pcd_index ++;
            string all_points_dir(string(string(ROOT_DIR) + "PCD/scans_") + to_string(pcd_index) + string(".pcd"));
            pcl::PCDWriter pcd_writer;
            cout << "current scan saved to /PCD/" << all_points_dir << endl;
            pcd_writer.writeBinary(all_points_dir, *pcl_wait_save);
            pcl_wait_save->clear();
            scan_wait_num = 0;
        }
    }
}

void publish_frame_body(const ros::Publisher & pubLaserCloudFull_body)
{
    int size = feats_undistort->points.size();
    PointCloudXYZI::Ptr laserCloudIMUBody(new PointCloudXYZI(size, 1));

    for (int i = 0; i < size; i++)
    {
        pointBodyLidarToIMU(&feats_undistort->points[i], \
                            &laserCloudIMUBody->points[i]);
    }

    sensor_msgs::PointCloud2 laserCloudmsg;
    pcl::toROSMsg(*laserCloudIMUBody, laserCloudmsg);
    laserCloudmsg.header.stamp = ros::Time().fromSec(lidar_end_time);
    laserCloudmsg.header.frame_id = "body";
    pubLaserCloudFull_body.publish(laserCloudmsg);
    publish_count -= PUBFRAME_PERIOD;
}

template<typename T>
void set_posestamp(T & out)
{
    if (!use_imu_as_input)
    {
        out.position.x = kf_output.x_.pos(0);
        out.position.y = kf_output.x_.pos(1);
        out.position.z = kf_output.x_.pos(2);
        out.orientation.x = kf_output.x_.rot.coeffs()[0];
        out.orientation.y = kf_output.x_.rot.coeffs()[1];
        out.orientation.z = kf_output.x_.rot.coeffs()[2];
        out.orientation.w = kf_output.x_.rot.coeffs()[3];
    }
    else
    {
        out.position.x = kf_input.x_.pos(0);
        out.position.y = kf_input.x_.pos(1);
        out.position.z = kf_input.x_.pos(2);
        out.orientation.x = kf_input.x_.rot.coeffs()[0];
        out.orientation.y = kf_input.x_.rot.coeffs()[1];
        out.orientation.z = kf_input.x_.rot.coeffs()[2];
        out.orientation.w = kf_input.x_.rot.coeffs()[3];
    }
}

void publish_odometry(const ros::Publisher & pubOdomAftMapped)
{
    odomAftMapped.header.frame_id = "camera_init";
    odomAftMapped.child_frame_id = "aft_mapped";
    if (publish_odometry_without_downsample)
    {
        odomAftMapped.header.stamp = ros::Time().fromSec(time_current);
    }
    else
    {
        odomAftMapped.header.stamp = ros::Time().fromSec(lidar_end_time);
    }
    set_posestamp(odomAftMapped.pose.pose);
    
    pubOdomAftMapped.publish(odomAftMapped);

    static tf::TransformBroadcaster br;
    tf::Transform                   transform;
    tf::Quaternion                  q;
    transform.setOrigin(tf::Vector3(odomAftMapped.pose.pose.position.x, \
                                    odomAftMapped.pose.pose.position.y, \
                                    odomAftMapped.pose.pose.position.z));
    q.setW(odomAftMapped.pose.pose.orientation.w);
    q.setX(odomAftMapped.pose.pose.orientation.x);
    q.setY(odomAftMapped.pose.pose.orientation.y);
    q.setZ(odomAftMapped.pose.pose.orientation.z);
    transform.setRotation( q );
    br.sendTransform( tf::StampedTransform( transform, odomAftMapped.header.stamp, "camera_init", "aft_mapped" ) );
}

void publish_path(const ros::Publisher pubPath)
{
    set_posestamp(msg_body_pose.pose);
    // msg_body_pose.header.stamp = ros::Time::now();
    msg_body_pose.header.stamp = ros::Time().fromSec(lidar_end_time);
    msg_body_pose.header.frame_id = "camera_init";
    static int jjj = 0;
    jjj++;
    // if (jjj % 2 == 0) // if path is too large, the rvis will crash
    {
        path.poses.emplace_back(msg_body_pose);
        pubPath.publish(path);
    }
}        

int main(int argc, char** argv)
{
    ros::init(argc, argv, "laserMapping");
    ros::NodeHandle nh("~");
    readParameters(nh);
    cout<<"lidar_type: "<<lidar_type<<endl;
    
    path.header.stamp    = ros::Time().fromSec(lidar_end_time);
    path.header.frame_id ="camera_init";

    /*** variables definition for counting ***/
    int frame_num = 0;
    double aver_time_consu = 0, aver_time_icp = 0, aver_time_match = 0, aver_time_incre = 0, aver_time_solve = 0, aver_time_propag = 0;
    
    /*** initialize variables ***/
    double FOV_DEG = (fov_deg + 10.0) > 179.9 ? 179.9 : (fov_deg + 10.0);
    double HALF_FOV_COS = cos((FOV_DEG) * 0.5 * PI_M / 180.0);

    memset(point_selected_surf, true, sizeof(point_selected_surf));
    downSizeFilterSurf.setLeafSize(filter_size_surf_min, filter_size_surf_min, filter_size_surf_min);
    downSizeFilterMap.setLeafSize(filter_size_map_min, filter_size_map_min, filter_size_map_min);
    if (!p_imu->UseLIInit)
    {
        Lidar_T_wrt_IMU<<VEC_FROM_ARRAY(extrinT);
        Lidar_R_wrt_IMU<<MAT_FROM_ARRAY(extrinR);
    }

    MatrixXd Jaco_rot(30000, 3);
    if (p_imu->UseLIInit)
    {
        p_imu->set_gyr_cov(V3D(li_init_gyr_cov, li_init_gyr_cov, li_init_gyr_cov));
        p_imu->set_acc_cov(V3D(li_init_acc_cov, li_init_acc_cov, li_init_acc_cov));

        Jaco_rot.setZero();
    }
    // if (extrinsic_est_en)
    // {
    //     if (!use_imu_as_input)
    //     {
    //         kf_output.x_.offset_R_L_I = Lidar_R_wrt_IMU;
    //         kf_output.x_.offset_T_L_I = Lidar_T_wrt_IMU;
    //     }
    //     else
    //     {
    //         kf_input.x_.offset_R_L_I = Lidar_R_wrt_IMU;
    //         kf_input.x_.offset_T_L_I = Lidar_T_wrt_IMU;
    //     }
    // }
    p_imu->lidar_type = p_pre->lidar_type = lidar_type;
    p_imu->imu_en = imu_en;
    if (GNSS_ENABLE)
    {
        std::copy(default_gnss_iono_params.begin(), default_gnss_iono_params.end(), 
            std::back_inserter(p_gnss->p_assign->latest_gnss_iono_params));
        p_gnss->Tex_imu_r << VEC_FROM_ARRAY(extrinT_gnss);
        p_gnss->Rex_imu_r << MAT_FROM_ARRAY(extrinR_gnss);
        p_gnss->gnss_ready = false; // gnss_quick_init; // edit
        p_gnss->nolidar = nolidar; // edit
        p_gnss->pre_integration->setnoise();

        if (p_gnss->p_assign->ephem_from_rinex)
        {
            p_gnss->p_assign->Ephemfromrinex(string("/home/joannahe-ssd/rosbag/gnss-lio/BRDM00DLR_S_20221870000_01D_MN.rnx"));
        }
    }

    kf_input.init_dyn_share_modified_2h(get_f_input, df_dx_input, h_model_input, h_model_GNSS_input);
    kf_output.init_dyn_share_modified_3h(get_f_output, df_dx_output, h_model_output, h_model_IMU_output, h_model_GNSS_output);
    Eigen::Matrix<double, 18, 18> P_init = MD(18, 18)::Identity() * 0.01;
    P_init.block<3, 3>(15, 15) = MD(3,3)::Identity() * 0.0001;
    P_init.block<6, 6>(9, 9) = MD(6,6)::Identity() * 0.001;
    // P_init.block<6, 6>(6, 6) = MD(6,6)::Identity() * 0.0001;
    kf_input.change_P(P_init);
    Eigen::Matrix<double, 24, 24> P_init_output = MD(24, 24)::Identity() * 0.01;
    P_init_output.block<3, 3>(15, 15) = MD(3,3)::Identity() * 0.0001;
    // P_init_output.block<6, 6>(6, 6) = MD(6,6)::Identity() * 0.0001;
    P_init_output.block<6, 6>(18, 18) = MD(6,6)::Identity() * 0.001;
    kf_input.change_P(P_init);
    kf_output.change_P(P_init_output);
    Eigen::Matrix<double, 18, 18> Q_input = process_noise_cov_input();
    Eigen::Matrix<double, 24, 24> Q_output = process_noise_cov_output();
    /*** debug record ***/
    FILE *fp;
    string pos_log_dir = root_dir + "/Log/pos_log.txt";
    fp = fopen(pos_log_dir.c_str(),"w");
    open_file();

    /*** ROS subscribe initialization ***/
    ros::Subscriber sub_pcl = p_pre->lidar_type == AVIA ? \
        nh.subscribe(lid_topic, 200000, livox_pcl_cbk) : \
        nh.subscribe(lid_topic, 200000, standard_pcl_cbk);
    ros::Subscriber sub_imu = nh.subscribe(imu_topic, 200000, imu_cbk);

    ros::Subscriber sub_ephem, sub_glo_ephem, sub_gnss_meas, sub_gnss_iono_params;
    ros::Subscriber sub_gnss_time_pluse_info, sub_local_trigger_info;
    ros::Subscriber sub_rtk_pvt_info, sub_rtk_lla_info;
    if (GNSS_ENABLE)
    {
        sub_ephem = nh.subscribe(gnss_ephem_topic, 1000, gnss_ephem_callback);
        sub_glo_ephem = nh.subscribe(gnss_glo_ephem_topic, 1000, gnss_glo_ephem_callback);
        sub_gnss_meas = nh.subscribe(gnss_meas_topic, 1000, gnss_meas_callback);
        sub_gnss_iono_params = nh.subscribe(gnss_iono_params_topic, 1000, gnss_iono_params_callback);
        sub_rtk_pvt_info = nh.subscribe(rtk_pvt_topic, 1000, rtk_pvt_callback);
        sub_rtk_lla_info = nh.subscribe(rtk_lla_topic, 1000, rtk_lla_callback);

        if (gnss_local_online_sync)
        {
            sub_gnss_time_pluse_info = nh.subscribe(gnss_tp_info_topic, 100, 
                gnss_tp_info_callback);
            sub_local_trigger_info = nh.subscribe(local_trigger_info_topic, 100, 
                local_trigger_info_callback);
        }
        else
        {
            time_diff_gnss_local = gnss_local_time_diff; // 18.0
            p_gnss->inputGNSSTimeDiff(time_diff_gnss_local);
            time_diff_valid = true;
        }
    }
    else
    {
        sub_rtk_pvt_info = nh.subscribe(rtk_pvt_topic, 100, rtk_pvt_callback);
    }

    ros::Publisher pubLaserCloudFullRes = nh.advertise<sensor_msgs::PointCloud2>
            ("/cloud_registered", 1000);
    ros::Publisher pubLaserCloudFullRes_body = nh.advertise<sensor_msgs::PointCloud2>
            ("/cloud_registered_body", 1000);
    ros::Publisher pubLaserCloudEffect  = nh.advertise<sensor_msgs::PointCloud2>
            ("/cloud_effected", 1000);
    ros::Publisher pubLaserCloudMap = nh.advertise<sensor_msgs::PointCloud2>
            ("/Laser_map", 1000);
    ros::Publisher pubOdomAftMapped = nh.advertise<nav_msgs::Odometry> 
            ("/aft_mapped_to_init", 1000);
    ros::Publisher pubPath          = nh.advertise<nav_msgs::Path> 
            ("/path", 1000);
    ros::Publisher plane_pub = nh.advertise<visualization_msgs::Marker>
            ("/planner_normal", 1000);
//------------------------------------------------------------------------------------------------------
    signal(SIGINT, SigHandle);
    ros::Rate rate(5000);
    bool status = ros::ok();
    while (status)
    {
        if (flg_exit) break;
        ros::spinOnce();
        if(sync_packages(Measures, p_gnss->gnss_msg)) 
        {
            if (flg_reset)
            {
                ROS_WARN("reset when rosbag play back");
                p_imu->Reset();
                feats_undistort.reset(new PointCloudXYZI());
                if (use_imu_as_input)
                {
                    state_in = kf_input.get_x();
                    state_in = state_input();
                    kf_input.change_P(P_init);
                }
                else
                {
                    state_out = kf_output.get_x();
                    state_out = state_output();
                    kf_output.change_P(P_init_output);
                }
                is_first_gnss = true;
                flg_first_scan = true;
                flg_reset = false;
                init_map = false;
                if (p_imu->UseLIInit)
                {
                    cut_frame_init = true; data_accum_finished = false; data_accum_start = false; online_calib_finish = false; refine_print = false;
                    frame_num_init = 0;
                    move_start_time = 0.0; online_calib_starts_time = 0.0;
                    Jaco_rot.setZero();
                }
                
                {
                    // ikdtree = KD_TREE<PointType>();
                    delete ikdtree.Root_Node;
                    ikdtree.Root_Node = nullptr;                
                }
                flg_reset = false;
            }

            if (flg_first_scan)
            {
                first_lidar_time = Measures.lidar_beg_time;
                flg_first_scan = false;
                if (first_imu_time < 1)
                {
                    first_imu_time = imu_next.header.stamp.toSec();
                    printf("first imu time: %f\n", first_imu_time);
                }
                time_current = 0.0;
                if(imu_en)
                {
                    // imu_next = *(imu_deque.front());
                    kf_input.x_.gravity << VEC_FROM_ARRAY(gravity);
                    kf_output.x_.gravity << VEC_FROM_ARRAY(gravity);
                    // kf_output.x_.acc << VEC_FROM_ARRAY(gravity);
                    // kf_output.x_.acc *= -1; 

                    if (!nolidar)
                    {
                        while (Measures.lidar_beg_time > imu_next.header.stamp.toSec()) // if it is needed for the new map?
                        {
                            imu_deque.pop_front();
                            if (imu_deque.empty())
                            {
                                break;
                            }
                            imu_last = imu_next;
                            imu_next = *(imu_deque.front());
                            // imu_deque.pop();
                        }
                    }
                }
                else
                {
                    kf_input.x_.gravity << VEC_FROM_ARRAY(gravity_init);
                    kf_output.x_.gravity << VEC_FROM_ARRAY(gravity_init);
                    kf_output.x_.acc << VEC_FROM_ARRAY(gravity_init);
                    kf_output.x_.acc *= -1; 
                    p_imu->imu_need_init_ = false;
                    p_imu->after_imu_init_ = true;
                }  
                if (!p_gnss->gnss_online_init && GNSS_ENABLE)
                {   
                    // p_gnss->gnss_ready = true;
                    // p_gnss->gtSAMgraphMade = true;
                    set_gnss_offline_init(false);
                }         
            }

            double t0,t1,t2,t3,t4,t5,match_start, solve_start;
            match_time = 0;
            solve_time = 0;
            propag_time = 0;
            update_time = 0;
            t0 = omp_get_wtime();
            
            /*** Segment the map in lidar FOV ***/
            lasermap_fov_segment();
            /*** downsample the feature points in a scan ***/
            t1 = omp_get_wtime();
            p_imu->Process(Measures, feats_undistort);
            if(space_down_sample)
            {
                downSizeFilterSurf.setInputCloud(feats_undistort);
                downSizeFilterSurf.filter(*feats_down_body);
                sort(feats_down_body->points.begin(), feats_down_body->points.end(), time_list); 
            }
            else
            {
                feats_down_body = Measures.lidar;
                sort(feats_down_body->points.begin(), feats_down_body->points.end(), time_list); 
            }
            if (!nolidar)
            {
                time_seq = time_compressing<int>(feats_down_body);
                feats_down_size = feats_down_body->points.size();
            }
            else
            {
                time_seq.clear();
            }
         
            if (!p_imu->UseLIInit && !p_imu->after_imu_init_)
            {
                if (!p_imu->imu_need_init_)
                { 
                    V3D tmp_gravity = - p_imu->mean_acc / acc_norm * G_m_s2;
                    // V3D tmp_gravity << VEC_FROM_ARRAY(gravity_init);
                    M3D rot_init;
                    p_imu->Set_init(tmp_gravity, rot_init);  
                    kf_input.x_.rot = kf_output.x_.rot = rot_init;
                    kf_input.x_.rot.normalize();
                    kf_output.x_.rot.normalize();
                    kf_output.x_.acc = - rot_init.transpose() * kf_output.x_.gravity;
                }
                continue;
            }
            /*** initialize the map kdtree ***/
            if(!init_map && !nolidar && !lose_lid)
            {
                if(ikdtree.Root_Node == nullptr) //
                // if(feats_down_size > 5)
                {
                    ikdtree.set_downsample_param(filter_size_map_min);
                }
                feats_down_world->resize(feats_undistort->size());
                for(int i = 0; i < feats_undistort->size(); i++)
                {
                    if (p_imu->UseLIInit)
                    {
                        p_imu->pointBodyToWorld_li_init(&(feats_undistort->points[i]), &(feats_down_world->points[i]));
                    }
                    else
                    {
                        pointBodyToWorld(&(feats_undistort->points[i]), &(feats_down_world->points[i]));
                    }
                }

                for (size_t i = 0; i < feats_down_world->size(); i++) 
                {
                    init_feats_world->points.emplace_back(feats_down_world->points[i]);
                }
                if(init_feats_world->size() < init_map_size) 
                {init_map = false;}
                else
                {
                    ikdtree.Build(init_feats_world->points); 
                    init_map = true;
                    publish_init_kdtree(pubLaserCloudMap); //(pubLaserCloudFullRes);
                }
                continue;
            }

            /*** ICP and Kalman filter update ***/
            normvec->resize(feats_down_size);
            feats_down_world->resize(feats_down_size);

            Nearest_Points.resize(feats_down_size);
            if (p_imu->UseLIInit)
            {
                LI_Init_update();
                /*** add the feature points to map kdtree ***/
                map_incremental();

                // if (scan_pub_en || pcd_save_en)      publish_frame_world(pubLaserCloudFullRes);
                // if (scan_pub_en && scan_body_pub_en) publish_frame_body(pubLaserCloudFullRes_body);
                frame_num_init ++;
            }

            if (!p_imu->LI_init_done && !data_accum_start && p_imu->state_LI_Init.pos.norm() > 0.05) {
                printf(BOLDCYAN "[Initialization] Movement detected, data accumulation starts.\n\n\n\n\n" RESET);
                data_accum_start = true;
                move_start_time = lidar_end_time;
            }
 
            if (!p_imu->LI_init_done && !data_accum_finished && data_accum_start) {
                //Push Lidar's Angular velocity and linear velocity
                Init_LI->push_Lidar_CalibState(p_imu->state_LI_Init.rot.normalized().toRotationMatrix(), p_imu->state_LI_Init.bg, p_imu->state_LI_Init.vel, lidar_end_time);
                //Data Accumulation Sufficience Appraisal
                data_accum_finished = Init_LI->data_sufficiency_assess(Jaco_rot, frame_num_init, p_imu->state_LI_Init.bg,
                                                                       orig_odom_freq, cut_frame_num);

                if (data_accum_finished) {
                    Init_LI->LI_Initialization(orig_odom_freq, cut_frame_num, timediff_imu_wrt_lidar, move_start_time);

                    online_calib_starts_time = lidar_end_time;
                    LI_Init_set();
                }
            }

            t2 = omp_get_wtime();
            
            if (p_imu->UseLIInit && p_imu->imu_need_init_ && !lose_lid)
            {
                continue;
            }
            /*** iterated state estimation ***/
            crossmat_list.reserve(feats_down_size);
            pbody_list.reserve(feats_down_size);
            // pbody_ext_list.reserve(feats_down_size);
                          
            for (size_t i = 0; i < feats_down_body->size(); i++)
            {
                V3D point_this(feats_down_body->points[i].x,
                            feats_down_body->points[i].y,
                            feats_down_body->points[i].z);
                pbody_list[i]=point_this;
                // if (extrinsic_est_en)
                // {
                //     if (!use_imu_as_input)
                //     {
                //         point_this = kf_output.x_.offset_R_L_I * point_this + kf_output.x_.offset_T_L_I;
                //     }
                //     else
                //     {
                //         point_this = kf_input.x_.offset_R_L_I * point_this + kf_input.x_.offset_T_L_I;
                //     }
                // }
                // else
                {
                    point_this = Lidar_R_wrt_IMU * point_this + Lidar_T_wrt_IMU;
                }
                M3D point_crossmat;
                point_crossmat << SKEW_SYM_MATRX(point_this);
                crossmat_list[i]=point_crossmat;
            }
            if (!use_imu_as_input)
            {     
                bool imu_upda_cov = false;
                effct_feat_num = 0;
                /**** point by point update ****/
                if (time_seq.size() > 0 || !GNSS_ENABLE)
                {
                double pcl_beg_time = Measures.lidar_beg_time;
                idx = -1;
                for (k = 0; k < time_seq.size(); k++)
                {
                    PointType &point_body  = feats_down_body->points[idx+time_seq[k]];

                    time_current = point_body.curvature / 1000.0 + pcl_beg_time;

                    if (is_first_frame)
                    {
                        if(imu_en)
                        {
                            while (time_current > imu_next.header.stamp.toSec())
                            {
                                imu_deque.pop_front();
                                if (imu_deque.empty()) break;
                                imu_last = imu_next;
                                imu_next = *(imu_deque.front());
                            }
                            angvel_avr<<imu_last.angular_velocity.x, imu_last.angular_velocity.y, imu_last.angular_velocity.z;
                            acc_avr   <<imu_last.linear_acceleration.x, imu_last.linear_acceleration.y, imu_last.linear_acceleration.z;
                        }
                        if (GNSS_ENABLE)
                        {
                            std::vector<Eigen::Vector3d>().swap(p_gnss->norm_vec_holder);
                            acc_avr_norm = acc_avr * G_m_s2 / acc_norm;
                            p_gnss->pre_integration->repropagate(kf_output.x_.ba, kf_output.x_.bg);
                            p_gnss->pre_integration->setacc0gyr0(acc_avr_norm, angvel_avr);
                        }
                        is_first_frame = false;
                        imu_upda_cov = true;
                        time_update_last = time_current;
                        time_predict_last_const = time_current;
                    }
                    if(imu_en && !imu_deque.empty())
                    {
                        bool last_imu = imu_next.header.stamp.toSec() == imu_deque.front()->header.stamp.toSec();
                        while (imu_next.header.stamp.toSec() < time_predict_last_const && !imu_deque.empty())
                        {
                            if (!last_imu)
                            {
                                imu_last = imu_next;
                                imu_next = *(imu_deque.front());
                                break;
                            }
                            else
                            {
                                imu_deque.pop_front();
                                if (imu_deque.empty()) break;
                                imu_last = imu_next;
                                imu_next = *(imu_deque.front());
                            }
                        }
                        bool imu_comes = time_current > imu_next.header.stamp.toSec();
                        while (imu_comes) 
                        {
                            if (!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                                while (time2sec(gnss_cur[0]->time) - gnss_local_time_diff < time_predict_last_const)
                                {
                                    p_gnss->gnss_msg.pop();
                                    if(!p_gnss->gnss_msg.empty())
                                    {
                                        gnss_cur = p_gnss->gnss_msg.front();
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                                while ((imu_next.header.stamp.toSec() > time2sec(gnss_cur[0]->time) - gnss_local_time_diff) && (time2sec(gnss_cur[0]->time) - gnss_local_time_diff >= time_predict_last_const))
                                {
                                    double dt = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_predict_last_const;
                                    double dt_cov = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_update_last;

                                    if (p_gnss->gnss_ready)
                                    {
                                        if (dt_cov > 0.0)
                                        {
                                            kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                        }
                                        kf_output.predict(dt, Q_output, input_in, true, false);
                                        p_gnss->pre_integration->push_back(dt, kf_output.x_.acc + kf_output.x_.ba, kf_output.x_.omg + kf_output.x_.bg); // acc_avr, angvel_avr); 
                                        p_gnss->processIMUOutput(dt, kf_output.x_.acc, kf_output.x_.omg);
                                        time_predict_last_const = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                        time_update_last = time_predict_last_const;
                                        p_gnss->processGNSS(gnss_cur, kf_output.x_);
                                        update_gnss = p_gnss->Evaluate(kf_output.x_);
                                        if (!p_gnss->gnss_ready)
                                        {
                                            flg_reset = true;
                                            p_gnss->gnss_msg.pop();
                                            if(!p_gnss->gnss_msg.empty())
                                            {
                                                gnss_cur = p_gnss->gnss_msg.front();
                                            }
                                            break; // ?
                                        }

                                        if (update_gnss)
                                        {
                                            kf_output.update_iterated_dyn_share_GNSS();
                                            cout_state_to_file();
                                        }
                                    }
                                    else
                                    {
                                        if (dt_cov > 0.0)
                                        {
                                            kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                        }
                                        
                                        kf_output.predict(dt_cov, Q_output, input_in, true, false);

                                        time_predict_last_const = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                        time_update_last = time_predict_last_const;
                                        state_out = kf_output.x_;
                                        state_out.rot = Rot_gnss_init * state_out.rot.normalized().toRotationMatrix();
                                        state_out.rot.normalize();
                                        state_out.pos = Rot_gnss_init * state_out.pos;
                                        state_out.vel = Rot_gnss_init * state_out.vel;
                                        p_gnss->processGNSS(gnss_cur, state_out);
                                    }
                                    p_gnss->gnss_msg.pop();
                                    if(!p_gnss->gnss_msg.empty())
                                    {
                                        gnss_cur = p_gnss->gnss_msg.front();
                                    }
                                    else
                                    {
                                        break;
                                    }
                                }
                            }
                            if (flg_reset)
                            {
                                break;
                            }
                            imu_upda_cov = true;
                            angvel_avr<<imu_next.angular_velocity.x, imu_next.angular_velocity.y, imu_next.angular_velocity.z;
                            acc_avr   <<imu_next.linear_acceleration.x, imu_next.linear_acceleration.y, imu_next.linear_acceleration.z;

                            /*** covariance update ***/
                            double dt = imu_next.header.stamp.toSec() - time_predict_last_const;
                            kf_output.predict(dt, Q_output, input_in, true, false);
                            time_predict_last_const = imu_next.header.stamp.toSec(); // big problem
                            // if (!imu_comes)
                            if (GNSS_ENABLE)
                            {
                                p_gnss->pre_integration->push_back(dt, kf_output.x_.acc + kf_output.x_.ba, kf_output.x_.omg + kf_output.x_.bg); // acc_avr, angvel_avr); 
                                p_gnss->processIMUOutput(dt, kf_output.x_.acc, kf_output.x_.omg);
                            }
                            {
                                double dt_cov = imu_next.header.stamp.toSec() - time_update_last; 

                                if (dt_cov > 0.0)
                                {
                                    time_update_last = imu_next.header.stamp.toSec();
                                    double propag_imu_start = omp_get_wtime();

                                    kf_output.predict(dt_cov, Q_output, input_in, false, true);

                                    propag_time += omp_get_wtime() - propag_imu_start;
                                    double solve_imu_start = omp_get_wtime();
                                    kf_output.update_iterated_dyn_share_IMU();
                                    solve_time += omp_get_wtime() - solve_imu_start;
                                }
                            }
                            imu_deque.pop_front();
                            if (imu_deque.empty()) break;
                            imu_last = imu_next;
                            imu_next = *(imu_deque.front());
                            imu_comes = time_current > imu_next.header.stamp.toSec();
                        }
                    }
                    if (flg_reset)
                    {
                        break;
                    }

                    if (!p_gnss->gnss_msg.empty())
                    {
                        gnss_cur = p_gnss->gnss_msg.front();
                        while ( time2sec(gnss_cur[0]->time) - gnss_local_time_diff < time_predict_last_const)
                        {
                            p_gnss->gnss_msg.pop();
                            if(!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                            }
                            else
                            {
                                break;
                            }
                        }
                        while (time_current > time2sec(gnss_cur[0]->time) - gnss_local_time_diff && time2sec(gnss_cur[0]->time) - gnss_local_time_diff >= time_predict_last_const)
                        {
                            double dt = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_predict_last_const;
                            double dt_cov = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_update_last;
                            // cout << "check gnss ready:" << p_gnss->gnss_ready << endl;
                            if (p_gnss->gnss_ready)
                            {
                                if (dt_cov > 0.0)
                                {
                                    kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                }
                                kf_output.predict(dt, Q_output, input_in, true, false);

                                p_gnss->pre_integration->push_back(dt, kf_output.x_.acc + kf_output.x_.ba, kf_output.x_.omg + kf_output.x_.bg); // acc_avr, angvel_avr); 
                                p_gnss->processIMUOutput(dt, kf_output.x_.acc, kf_output.x_.omg);

                                time_predict_last_const = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                time_update_last = time_predict_last_const;
                                p_gnss->processGNSS(gnss_cur, kf_output.x_);
                                update_gnss = p_gnss->Evaluate(kf_output.x_);
                                if (!p_gnss->gnss_ready)
                                {
                                    flg_reset = true;
                                    p_gnss->gnss_msg.pop();
                                    if(!p_gnss->gnss_msg.empty())
                                    {
                                        gnss_cur = p_gnss->gnss_msg.front();
                                    }
                                    break; // ?
                                }

                                if (update_gnss)
                                {
                                    kf_output.update_iterated_dyn_share_GNSS();
                                    cout_state_to_file();
                                }
                            }
                            else
                            {
                                if (dt_cov > 0.0)
                                {
                                    kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                }
                                kf_output.predict(dt, Q_output, input_in, false, true);
                                time_predict_last_const = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                time_update_last = time_predict_last_const;
                                state_out = kf_output.x_;
                                state_out.rot = Rot_gnss_init * state_out.rot.normalized().toRotationMatrix();
                                state_out.rot.normalize();
                                state_out.pos = Rot_gnss_init * state_out.pos;
                                state_out.vel = Rot_gnss_init * state_out.vel;
                                p_gnss->processGNSS(gnss_cur, state_out);
                            }
                            p_gnss->gnss_msg.pop();
                            if(!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                            }
                            else
                            {
                                break;
                            }
                        }
                    }
                    if (flg_reset)
                    {
                        break;
                    }
                    double dt = time_current - time_predict_last_const;
                    double propag_state_start = omp_get_wtime();
                    if(!prop_at_freq_of_imu)
                    {
                        double dt_cov = time_current - time_update_last;
                        if (dt_cov > 0.0)
                        {
                            kf_output.predict(dt_cov, Q_output, input_in, false, true);
                            time_update_last = time_current;   
                        }
                    }
                    kf_output.predict(dt, Q_output, input_in, true, false);
                    propag_time += omp_get_wtime() - propag_state_start;
                    time_predict_last_const = time_current;
                    if (GNSS_ENABLE)
                    {
                        p_gnss->pre_integration->push_back(dt, kf_output.x_.acc + kf_output.x_.ba, kf_output.x_.omg + kf_output.x_.bg); // acc_avr, angvel_avr); 
                        p_gnss->processIMUOutput(dt, kf_output.x_.acc, kf_output.x_.omg);
                    }
                    // if(k == 0)
                    // {
                    //     fout_imu_pbp << Measures.lidar_last_time - first_lidar_time << " " << imu_last.angular_velocity.x << " " << imu_last.angular_velocity.y << " " << imu_last.angular_velocity.z \
                    //             << " " << imu_last.linear_acceleration.x << " " << imu_last.linear_acceleration.y << " " << imu_last.linear_acceleration.z << endl;
                    // }
                    double t_update_start = omp_get_wtime();

                    if (feats_down_size < 1)
                    {
                        ROS_WARN("No point, skip this scan!\n");
                        idx += time_seq[k];
                        continue;
                    }
                    if (!kf_output.update_iterated_dyn_share_modified()) 
                    {
                        idx = idx+time_seq[k];
                        continue;
                    }
                    solve_start = omp_get_wtime();
                        
                    if (publish_odometry_without_downsample)
                    {
                        /******* Publish odometry *******/

                        publish_odometry(pubOdomAftMapped);
                        if (runtime_pos_log)
                        {
                            euler_cur = SO3ToEuler(kf_output.x_.rot);
                            fout_out << setw(20) << Measures.lidar_beg_time - first_lidar_time << " " << euler_cur.transpose() << " " << kf_output.x_.pos.transpose() << " " << kf_output.x_.vel.transpose() \
                            <<" "<<kf_output.x_.omg.transpose()<<" "<<kf_output.x_.acc.transpose()<<" "<<kf_output.x_.gravity.transpose()<<" "<<kf_output.x_.bg.transpose()<<" "<<kf_output.x_.ba.transpose()<<" "<<feats_undistort->points.size()<<endl;
                        }
                    }

                    for (int j = 0; j < time_seq[k]; j++)
                    {
                        PointType &point_body_j  = feats_down_body->points[idx+j+1];
                        PointType &point_world_j = feats_down_world->points[idx+j+1];
                        pointBodyToWorld(&point_body_j, &point_world_j);
                    }
                
                    solve_time += omp_get_wtime() - solve_start;
    
                    update_time += omp_get_wtime() - t_update_start;
                    idx += time_seq[k];
                    // cout << "pbp output effect feat num:" << effct_feat_num << endl;
                }
                }
                else
                {
                    p_gnss->nolidar_cur = true;
                    if (!imu_deque.empty())
                    { 
                        imu_last = imu_next;
                        imu_next = *(imu_deque.front());

                    while (imu_next.header.stamp.toSec() > time_current && ((imu_next.header.stamp.toSec() < imu_first_time + lidar_time_inte && nolidar) || (imu_next.header.stamp.toSec() < Measures.lidar_beg_time + lidar_time_inte && !nolidar)))
                    { // >= ?
                        if (is_first_frame)
                        {
                            if (!nolidar) std::vector<Eigen::Vector3d>().swap(p_gnss->norm_vec_holder);
                            if (!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                                double front_gnss_ts = time2sec(gnss_cur[0]->time); // take time
                                time_current = front_gnss_ts - gnss_local_time_diff;
                                while (imu_next.header.stamp.toSec() < time_current) // 0.05
                                {
                                    ROS_WARN("throw IMU, only should happen at the beginning 2510");
                                    imu_deque.pop_front();
                                    if (imu_deque.empty()) break;
                                    imu_last = imu_next;
                                    imu_next = *(imu_deque.front()); // could be used to initialize
                                }
                            }
                            else
                            {
                                if (nolidar)
                                {
                                    while (imu_next.header.stamp.toSec() < imu_first_time + lidar_time_inte)
                                    {
                                        // meas.imu.emplace_back(imu_deque.front()); should add to initialization
                                        imu_deque.pop_front();
                                        if(imu_deque.empty()) break;
                                        imu_last = imu_next;
                                        imu_next = *(imu_deque.front()); // could be used to initialize
                                    }
                                }
                                else
                                {
                                    while (imu_next.header.stamp.toSec() < Measures.lidar_beg_time + lidar_time_inte)
                                    {
                                        // meas.imu.emplace_back(imu_deque.front()); should add to initialization
                                        imu_deque.pop_front();
                                        if(imu_deque.empty()) break;
                                        imu_last = imu_next;
                                        imu_next = *(imu_deque.front());
                                    }
                                }
                                break;
                            }
                            angvel_avr<<imu_last.angular_velocity.x, imu_last.angular_velocity.y, imu_last.angular_velocity.z;
                                            
                            acc_avr   <<imu_last.linear_acceleration.x, imu_last.linear_acceleration.y, imu_last.linear_acceleration.z;

                            imu_upda_cov = true;
                            time_update_last = time_current;
                            time_predict_last_const = time_current;
                            p_gnss->pre_integration->repropagate(kf_output.x_.ba, kf_output.x_.bg);
                            acc_avr_norm = acc_avr * G_m_s2 / acc_norm;
                            p_gnss->pre_integration->setacc0gyr0(acc_avr_norm, angvel_avr); 

                            if (nolidar && !p_gnss->gnss_online_init)
                            {
                                if (!p_gnss->gnss_msg.empty())
                                {
                                    gnss_cur = p_gnss->gnss_msg.front();
                                    if (time2sec(gnss_cur[0]->time) < imu_first_time + lidar_time_inte + gnss_local_time_diff)
                                    {
                                        p_gnss->processGNSS(gnss_cur, kf_output.x_);
                                        if (1) // (p_gnss->gnss_meas_buf[0].size() > 4)
                                        {
                                            set_gnss_offline_init(true);
                                        }
                                    }
                                }
                            }
                            else
                            {
                                is_first_frame = false;
                            }
                        }
                        time_current = imu_next.header.stamp.toSec();

                        if (!is_first_frame)
                        {
                        if (!p_gnss->gnss_msg.empty())
                        {
                            gnss_cur = p_gnss->gnss_msg.front();
                            while ( time2sec(gnss_cur[0]->time) - gnss_local_time_diff < time_predict_last_const)
                            {
                                p_gnss->gnss_msg.pop();
                                if(!p_gnss->gnss_msg.empty())
                                {
                                    gnss_cur = p_gnss->gnss_msg.front();
                                }
                                else
                                {
                                    break;
                                }
                            }
                        while ((time_current > time2sec(gnss_cur[0]->time) - gnss_local_time_diff) && (time2sec(gnss_cur[0]->time) - gnss_local_time_diff >= time_predict_last_const))
                        {
                            double dt = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_predict_last_const;
                            double dt_cov = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_update_last;

                            if (p_gnss->gnss_ready)
                            {
                                if (dt_cov > 0.0)
                                {
                                    kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                    time_update_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                }
                                kf_output.predict(dt, Q_output, input_in, true, false);
                                p_gnss->pre_integration->push_back(dt, kf_output.x_.acc + kf_output.x_.ba, kf_output.x_.omg + kf_output.x_.bg); //acc_avr_norm, angvel_avr); 
                                // change to state_const.omg and state_const.acc? 

                                time_predict_last_const = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                p_gnss->processGNSS(gnss_cur, kf_output.x_);
                                update_gnss = p_gnss->Evaluate(kf_output.x_); 
                                if (!p_gnss->gnss_ready)
                                {
                                    flg_reset = true;
                                    p_gnss->gnss_msg.pop();
                                    if(!p_gnss->gnss_msg.empty())
                                    {
                                        gnss_cur = p_gnss->gnss_msg.front();
                                    }
                                    break; // ?
                                }
                                if (update_gnss)
                                {
                                    if (!nolidar)
                                    {
                                        kf_output.update_iterated_dyn_share_GNSS();
                                    }
                                    cout_state_to_file();
                                }
                            }
                            else
                            {
                                if (dt_cov > 0.0)
                                {
                                    kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                    time_update_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                }
                                kf_output.predict(dt, Q_output, input_in, true, false);
                                time_predict_last_const = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                p_gnss->processGNSS(gnss_cur, kf_output.x_);
                                if (p_gnss->gnss_ready)
                                {
                                    if (nolidar)
                                    {
                                        // Eigen::Matrix3d R_enu_local_;
                                        // R_enu_local_ = Eigen::AngleAxisd(p_gnss->yaw_enu_local, Eigen::Vector3d::UnitZ());
                                        kf_output.x_.pos = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(p_gnss->frame_num-1)).segment<3>(0); // p_gnss->anc_ecef - p_gnss->R_ecef_enu * R_enu_local_ * state_const.rot_end * p_gnss->Tex_imu_r;
                                        kf_output.x_.rot = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Rot3>(R(p_gnss->frame_num-1)).matrix(); // p_gnss->R_ecef_enu * R_enu_local_ * state_const.rot_end;
                                        kf_output.x_.rot.normalize();
                                        kf_output.x_.vel = p_gnss->p_assign->isamCurrentEstimate.at<gtsam::Vector12>(F(p_gnss->frame_num-1)).segment<3>(3); // p_gnss->R_ecef_enu * R_enu_local_ * state_const.vel_end; // Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                        kf_output.x_.ba = Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                        kf_output.x_.bg = Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                        kf_output.x_.omg = Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                        kf_output.x_.gravity = p_gnss->R_ecef_enu * kf_output.x_.gravity; // * R_enu_local_ 
                                        kf_output.x_.acc = kf_output.x_.rot.conjugate().normalized() * (-kf_output.x_.gravity); // R_ecef_enu * state.vel_end;
                                        
                                        kf_output.P_ = MD(24,24)::Identity() * INIT_COV;
                                    }
                                }
                            }
                            p_gnss->gnss_msg.pop();
                            if(!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                            }
                            else
                            {
                                break;
                            }
                        }
                        }
                        if (flg_reset)
                        {
                            break;
                        }
                        double dt = time_current - time_predict_last_const;
                        {
                            double dt_cov = time_current - time_update_last;
                            if (dt_cov > 0.0)
                            {
                                kf_output.predict(dt_cov, Q_output, input_in, false, true);
                                time_update_last = time_current;
                            }
                            kf_output.predict(dt, Q_output, input_in, true, false);

                            p_gnss->pre_integration->push_back(dt, kf_output.x_.acc + kf_output.x_.ba, kf_output.x_.omg + kf_output.x_.bg); // acc_avr_norm, angvel_avr); // 
                        }

                        time_predict_last_const = time_current;

                        angvel_avr<<imu_next.angular_velocity.x, imu_next.angular_velocity.y, imu_next.angular_velocity.z;
                        acc_avr   <<imu_next.linear_acceleration.x, imu_next.linear_acceleration.y, imu_next.linear_acceleration.z; 
                        acc_avr_norm = acc_avr * G_m_s2 / acc_norm;
                        kf_output.update_iterated_dyn_share_IMU();
                        imu_deque.pop_front();
                        if (imu_deque.empty()) break;
                        imu_last = imu_next;
                        imu_next = *(imu_deque.front());
                    }
                    else
                    {
                        imu_deque.pop_front();
                        if (imu_deque.empty()) break;
                        imu_last = imu_next;
                        imu_next = *(imu_deque.front());
                    }
                    }
                    }
                }
            }
            else
            {
                bool imu_prop_cov = false;
                effct_feat_num = 0;
                if (time_seq.size() > 0 || !GNSS_ENABLE)
                {
                double pcl_beg_time = Measures.lidar_beg_time;
                idx = -1;
                for (k = 0; k < time_seq.size(); k++)
                {
                    PointType &point_body  = feats_down_body->points[idx+time_seq[k]];
                    time_current = point_body.curvature / 1000.0 + pcl_beg_time;
                    if (is_first_frame)
                    {
                        while (time_current > imu_next.header.stamp.toSec()) 
                        {
                            imu_deque.pop_front();
                            if (imu_deque.empty()) break;
                            imu_last = imu_next;
                            imu_next = *(imu_deque.front());
                        }
                        imu_prop_cov = true;

                        is_first_frame = false;
                        t_last = time_current;
                        time_update_last = time_current; 
                        {
                            input_in.gyro<<imu_last.angular_velocity.x, imu_last.angular_velocity.y, imu_last.angular_velocity.z;                 
                            input_in.acc<<imu_last.linear_acceleration.x, imu_last.linear_acceleration.y, imu_last.linear_acceleration.z;
                            input_in.acc = input_in.acc * G_m_s2 / acc_norm;
                        }
                        if (GNSS_ENABLE)
                        {
                            p_gnss->pre_integration->repropagate(kf_input.x_.ba, kf_input.x_.bg);
                            p_gnss->pre_integration->setacc0gyr0(input_in.acc, input_in.gyro);
                            std::vector<Eigen::Vector3d>().swap(p_gnss->norm_vec_holder);
                        }
                    }
                    
                    while (time_current > imu_next.header.stamp.toSec()) // && !imu_deque.empty())
                    {
                        imu_deque.pop_front();
                        if (!p_gnss->gnss_msg.empty())
                        {
                            gnss_cur = p_gnss->gnss_msg.front();
                            while (time2sec(gnss_cur[0]->time) - gnss_local_time_diff < t_last)
                            {
                                p_gnss->gnss_msg.pop();
                                if(!p_gnss->gnss_msg.empty())
                                {
                                    gnss_cur = p_gnss->gnss_msg.front();
                                }
                                else
                                {
                                    break;
                                }
                            }
                            while ((imu_next.header.stamp.toSec() > time2sec(gnss_cur[0]->time) - gnss_local_time_diff) && (time2sec(gnss_cur[0]->time) - gnss_local_time_diff >= t_last))
                            {
                                double dt = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - t_last;
                                double dt_cov = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_update_last;
                                // cout << "check gnss ready:" << p_gnss->gnss_ready << endl;
                                if (p_gnss->gnss_ready)
                                {
                                    if (dt_cov > 0.0)
                                    {
                                        kf_input.predict(dt_cov, Q_input, input_in, false, true);
                                        time_update_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff; 
                                    }
                                    kf_input.predict(dt_cov, Q_input, input_in, false, true);

                                    p_gnss->pre_integration->push_back(dt, input_in.acc, input_in.gyro);
                                    p_gnss->processIMU(dt, input_in.acc, input_in.gyro);

                                    t_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                    p_gnss->processGNSS(gnss_cur, kf_input.x_, input_in.gyro);
                                    update_gnss = p_gnss->Evaluate(kf_input.x_, input_in.gyro);
                                    if (!p_gnss->gnss_ready)
                                    {
                                        flg_reset = true;
                                        p_gnss->gnss_msg.pop();
                                        if(!p_gnss->gnss_msg.empty())
                                        {
                                            gnss_cur = p_gnss->gnss_msg.front();
                                        }
                                        break; // ?
                                    }
                                    if (update_gnss)
                                    {
                                        kf_input.update_iterated_dyn_share_GNSS();
                                        cout_state_to_file();
                                    }
                                }
                                else
                                {
                                    if (dt_cov > 0.0)
                                    {
                                        kf_input.predict(dt_cov, Q_input, input_in, false, true);
                                        time_update_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff; //time_current;
                                    }
                                    kf_input.predict(dt, Q_input, input_in, true, false);
                                    t_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                    state_in = kf_input.x_;
                                    state_in.rot = Rot_gnss_init * state_in.rot.normalized().toRotationMatrix();
                                    state_in.rot.normalize();
                                    state_in.pos = Rot_gnss_init * state_in.pos;
                                    state_in.vel = Rot_gnss_init * state_in.vel;
                                    p_gnss->processGNSS(gnss_cur, state_in, input_in.gyro);
                                }
                                p_gnss->gnss_msg.pop();
                                if(!p_gnss->gnss_msg.empty())
                                {
                                    gnss_cur = p_gnss->gnss_msg.front();
                                }
                                else
                                {
                                    break;
                                }
                            }
                        }
                        if (flg_reset)
                        {
                            break;
                        }
                        
                        input_in.gyro<<imu_last.angular_velocity.x, imu_last.angular_velocity.y, imu_last.angular_velocity.z;
                        input_in.acc <<imu_last.linear_acceleration.x, imu_last.linear_acceleration.y, imu_last.linear_acceleration.z; 
                        input_in.acc    = input_in.acc * G_m_s2 / acc_norm; 
                        double dt = imu_last.header.stamp.toSec() - t_last;

                        double dt_cov = imu_last.header.stamp.toSec() - time_update_last;
                        if (p_gnss->gnss_ready)
                        {
                            p_gnss->pre_integration->push_back(dt, input_in.acc, input_in.gyro);
                            p_gnss->processIMU(dt, input_in.acc, input_in.gyro);
                        }
                        if (dt_cov > 0.0)
                        {
                            kf_input.predict(dt_cov, Q_input, input_in, false, true); 
                            time_update_last = imu_last.header.stamp.toSec(); //time_current;
                        }
                        kf_input.predict(dt, Q_input, input_in, true, false); 
                        t_last = imu_last.header.stamp.toSec();
                        imu_prop_cov = true;

                        if (imu_deque.empty()) break;
                        imu_last = imu_next;
                        imu_next = *(imu_deque.front());
                        // imu_upda_cov = true;
                    }     
                    if (flg_reset)
                    {
                        break;
                    }     
                    if (!p_gnss->gnss_msg.empty())
                    {
                        gnss_cur = p_gnss->gnss_msg.front();
                        while (time2sec(gnss_cur[0]->time) - gnss_local_time_diff < t_last)
                        {
                            p_gnss->gnss_msg.pop();
                            if(!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                            }
                            else
                            {
                                break;
                            }
                        }
                    while ((time_current > time2sec(gnss_cur[0]->time) - gnss_local_time_diff) && (time2sec(gnss_cur[0]->time) - gnss_local_time_diff >= t_last))
                    {
                        double dt = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - t_last;
                        double dt_cov = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_update_last;
 
                        if (p_gnss->gnss_ready)
                        {
                            if (dt_cov > 0.0)
                            {
                                kf_input.predict(dt_cov, Q_input, input_in, false, true); 
                                time_update_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff; //time_current;
                            }
                            kf_input.predict(dt, Q_input, input_in, true, false); 

                            p_gnss->pre_integration->push_back(dt, input_in.acc, input_in.gyro);
                            p_gnss->processIMU(dt, input_in.acc, input_in.gyro);
                            t_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                            p_gnss->processGNSS(gnss_cur, kf_input.x_, input_in.gyro);
                            update_gnss = p_gnss->Evaluate(kf_input.x_, input_in.gyro);
                            if (!p_gnss->gnss_ready)
                            {
                                flg_reset = true;
                                p_gnss->gnss_msg.pop();
                                if(!p_gnss->gnss_msg.empty())
                                {
                                    gnss_cur = p_gnss->gnss_msg.front();
                                }
                                break; // ?
                            }
                            if (update_gnss)
                            {
                                kf_input.update_iterated_dyn_share_GNSS();
                                cout_state_to_file();
                            }
                        }
                        else
                        {
                            if (dt_cov > 0.0)
                            {                                
                                kf_input.predict(dt_cov, Q_input, input_in, false, true); 
                                time_update_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff; //time_current;
                            }
                            kf_input.predict(dt, Q_input, input_in, true, false); 
                            t_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                            state_in = kf_input.x_;
                            state_in.rot = Rot_gnss_init * state_in.rot.normalized().toRotationMatrix();
                            state_in.rot.normalize();
                            state_in.pos = Rot_gnss_init * state_in.pos;
                            state_in.vel = Rot_gnss_init * state_in.vel;
                            p_gnss->processGNSS(gnss_cur, state_in, input_in.gyro);
                        }
                        p_gnss->gnss_msg.pop();
                        if(!p_gnss->gnss_msg.empty())
                        {
                            gnss_cur = p_gnss->gnss_msg.front();
                        }
                        else
                        {
                            break;
                        }
                    } 
                    }
                    if (flg_reset)
                    {
                        break;
                    }
                    double dt = time_current - t_last;
                    t_last = time_current;
                    double propag_start = omp_get_wtime();
                    
                    if(!prop_at_freq_of_imu)
                    {   
                        double dt_cov = time_current - time_update_last;
                        if (dt_cov > 0.0)
                        {    
                            kf_input.predict(dt_cov, Q_input, input_in, false, true); 
                            time_update_last = time_current; 
                        }
                    }
                    kf_input.predict(dt, Q_input, input_in, true, false); 

                    propag_time += omp_get_wtime() - propag_start;
                    if (GNSS_ENABLE)
                    {
                        p_gnss->pre_integration->push_back(dt, input_in.acc, input_in.gyro);
                        p_gnss->processIMU(dt, input_in.acc, input_in.gyro);
                    }
                    // if(k == 0)
                    // {
                    //     fout_imu_pbp << Measures.lidar_last_time - first_lidar_time << " " << imu_last.angular_velocity.x << " " << imu_last.angular_velocity.y << " " << imu_last.angular_velocity.z \
                    //             << " " << imu_last.linear_acceleration.x << " " << imu_last.linear_acceleration.y << " " << imu_last.linear_acceleration.z << endl;
                    // }

                    double t_update_start = omp_get_wtime();
                    
                    if (feats_down_size < 1)
                    {
                        ROS_WARN("No point, skip this scan!\n");

                        idx += time_seq[k];
                        continue;
                    }
                    if (!kf_input.update_iterated_dyn_share_modified()) 
                    {
                        idx = idx+time_seq[k];
                        continue;
                    }

                    solve_start = omp_get_wtime();

                    if (publish_odometry_without_downsample)
                    {
                        /******* Publish odometry *******/

                        publish_odometry(pubOdomAftMapped);
                        if (runtime_pos_log)
                        {
                            euler_cur = SO3ToEuler(kf_input.x_.rot);
                            fout_out << setw(20) << Measures.lidar_beg_time - first_lidar_time << " " << euler_cur.transpose() << " " << kf_input.x_.pos.transpose() << " " << kf_input.x_.vel.transpose() \
                            <<" "<<kf_input.x_.bg.transpose()<<" "<<kf_input.x_.ba.transpose()<<" "<<kf_input.x_.gravity.transpose()<<" "<<feats_undistort->points.size()<<endl;
                        }
                    }

                    for (int j = 0; j < time_seq[k]; j++)
                    {
                        PointType &point_body_j  = feats_down_body->points[idx+j+1];
                        PointType &point_world_j = feats_down_world->points[idx+j+1];
                        pointBodyToWorld(&point_body_j, &point_world_j); 
                    }
                    solve_time += omp_get_wtime() - solve_start;
                
                    update_time += omp_get_wtime() - t_update_start;
                    idx = idx + time_seq[k];
                }  
                }
                else
                {
                    p_gnss->nolidar_cur = true;
                    if (!imu_deque.empty())
                    { 
                    imu_last = imu_next;
                    imu_next = *(imu_deque.front());
                    while (imu_next.header.stamp.toSec() > time_current && ((imu_next.header.stamp.toSec() < imu_first_time + lidar_time_inte && nolidar )|| (imu_next.header.stamp.toSec() < Measures.lidar_beg_time + lidar_time_inte && !nolidar)))
                    { // >= ?
                        if (is_first_frame)
                        {
                            if (!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                                double front_gnss_ts = time2sec(gnss_cur[0]->time); // take time
                                time_current = front_gnss_ts - gnss_local_time_diff;
                                while (imu_next.header.stamp.toSec() < time_current) // 0.05
                                {
                                    ROS_WARN("throw IMU, only should happen at the beginning 2510");
                                    imu_deque.pop_front();
                                    if (imu_deque.empty()) break;
                                    imu_last = imu_next;
                                    imu_next = *(imu_deque.front()); // could be used to initialize
                                }
                            }
                            else
                            {
                                if (nolidar)
                                {
                                    while (imu_next.header.stamp.toSec() < imu_first_time + lidar_time_inte)
                                    {
                                        imu_deque.pop_front();
                                        if(imu_deque.empty()) break;
                                        imu_last = imu_next;
                                        imu_next = *(imu_deque.front()); // could be used to initialize
                                    }
                                }
                                else
                                {
                                    while (imu_next.header.stamp.toSec() < Measures.lidar_beg_time + lidar_time_inte)
                                    {
                                        imu_deque.pop_front();
                                        if(imu_deque.empty()) break;
                                        imu_last = imu_next;
                                        imu_next = *(imu_deque.front());
                                    }
                                }
                                
                                break;
                            }
                            imu_prop_cov = true;
                            
                            t_last = time_current;
                            time_update_last = time_current; 
                            input_in.gyro<<imu_last.angular_velocity.x, imu_last.angular_velocity.y, imu_last.angular_velocity.z;
                            input_in.acc   <<imu_last.linear_acceleration.x, imu_last.linear_acceleration.y, imu_last.linear_acceleration.z;
                            input_in.acc = input_in.acc * G_m_s2 / acc_norm;
                            if (!nolidar) std::vector<Eigen::Vector3d>().swap(p_gnss->norm_vec_holder);
                            p_gnss->pre_integration->repropagate(kf_input.x_.ba, kf_input.x_.bg);
                            p_gnss->pre_integration->setacc0gyr0(input_in.acc, input_in.gyro);
                            if (nolidar && !p_gnss->gnss_online_init) // no meaning
                            {
                                if (!p_gnss->gnss_msg.empty())
                                {
                                    gnss_cur = p_gnss->gnss_msg.front();
                                    if (time2sec(gnss_cur[0]->time) < imu_first_time + lidar_time_inte + gnss_local_time_diff)
                                    {
                                        p_gnss->processGNSS(gnss_cur, kf_input.x_, input_in.gyro);
                                        if (1) // (p_gnss->gnss_meas_buf[0].size() > 4)
                                        {
                                            set_gnss_offline_init(true);
                                        }
                                    }
                                }
                            }
                            else
                            {
                                is_first_frame = false;
                            }
                        }
                        time_current = imu_next.header.stamp.toSec();

                        if (!is_first_frame)
                        {
                        if (!p_gnss->gnss_msg.empty())
                        {
                            gnss_cur = p_gnss->gnss_msg.front();
                            while (time2sec(gnss_cur[0]->time) - gnss_local_time_diff < t_last)
                            {
                                p_gnss->gnss_msg.pop();
                                if(!p_gnss->gnss_msg.empty())
                                {
                                    gnss_cur = p_gnss->gnss_msg.front();
                                }
                                else
                                {
                                    break;
                                }
                            }
                        while ((time_current > time2sec(gnss_cur[0]->time) - gnss_local_time_diff) && (time2sec(gnss_cur[0]->time) - gnss_local_time_diff >= t_last))
                        {
                            double dt = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - t_last;
                            double dt_cov = time2sec(gnss_cur[0]->time) - gnss_local_time_diff - time_update_last;

                            if (p_gnss->gnss_ready)
                            {
                                if (dt_cov > 0.0)
                                {
                                    kf_input.predict(dt_cov, Q_input, input_in, false, true);
                                    time_update_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff; //time_current;
                                }
                                kf_input.predict(dt, Q_input, input_in, true, false);

                                p_gnss->pre_integration->push_back(dt, input_in.acc, input_in.gyro);

                                t_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                p_gnss->processGNSS(gnss_cur, kf_input.x_, input_in.gyro);
                                update_gnss = p_gnss->Evaluate(kf_input.x_, input_in.gyro);
                                if (!p_gnss->gnss_ready)
                                {
                                    flg_reset = true;
                                    p_gnss->gnss_msg.pop();
                                    if(!p_gnss->gnss_msg.empty())
                                    {
                                        gnss_cur = p_gnss->gnss_msg.front();
                                    }
                                    break; // ?
                                }
                                if (update_gnss)
                                {
                                    if (!nolidar)
                                    {
                                        kf_input.update_iterated_dyn_share_GNSS();
                                    }
                                    cout_state_to_file();
                                }
                            }
                            else
                            {
                                if (dt_cov > 0.0)
                                {
                                    kf_input.predict(dt_cov, Q_input, input_in, false, true);
                                    time_update_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff; //time_current;
                                }
                                kf_input.predict(dt, Q_input, input_in, true, false);
                                t_last = time2sec(gnss_cur[0]->time) - gnss_local_time_diff;
                                p_gnss->processGNSS(gnss_cur, kf_input.x_, input_in.gyro);
                                if (p_gnss->gnss_ready)
                                {
                                    if (nolidar)
                                    {
                                        Eigen::Matrix3d R_enu_local_;
                                        R_enu_local_ = Eigen::AngleAxisd(p_gnss->yaw_enu_local, Eigen::Vector3d::UnitZ());
                                        kf_input.x_.pos = p_gnss->anc_ecef - p_gnss->R_ecef_enu * R_enu_local_ * kf_input.x_.rot.normalized().toRotationMatrix() * p_gnss->Tex_imu_r;
                                        kf_input.x_.rot = p_gnss->R_ecef_enu * R_enu_local_ * kf_input.x_.rot.normalized().toRotationMatrix();
                                        kf_input.x_.rot.normalize();
                                        kf_input.x_.vel = p_gnss->R_ecef_enu * R_enu_local_ * kf_input.x_.vel; //Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                        kf_input.x_.ba = Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                        kf_input.x_.bg = Eigen::Vector3d::Zero(); // R_ecef_enu * state.vel_end;
                                        kf_input.x_.gravity = p_gnss->R_ecef_enu * kf_input.x_.gravity; // * R_enu_local_
                                        kf_input.P_ = MD(18,18)::Identity() * INIT_COV;
                                    }
                                }
                            }
                            p_gnss->gnss_msg.pop();
                            if(!p_gnss->gnss_msg.empty())
                            {
                                gnss_cur = p_gnss->gnss_msg.front();
                            }
                            else
                            {
                                break;
                            }
                        }
                        }
                        if (flg_reset)
                        {
                            break;
                        }
                        double dt = time_current - t_last;

                        double dt_cov = time_current - time_update_last;
                        if (dt_cov > 0.0)
                        {        
                            kf_input.predict(dt_cov, Q_input, input_in, false, true);
                            time_update_last = imu_next.header.stamp.toSec(); //time_current;
                        }
                        kf_input.predict(dt, Q_input, input_in, true, false);

                        t_last = imu_next.header.stamp.toSec();
                    
                        if (p_gnss->gnss_ready)
                        {
                            p_gnss->pre_integration->push_back(dt, input_in.acc, input_in.gyro);
                        }
                        input_in.gyro<<imu_next.angular_velocity.x, imu_next.angular_velocity.y, imu_next.angular_velocity.z;
                        input_in.acc<<imu_next.linear_acceleration.x, imu_next.linear_acceleration.y, imu_next.linear_acceleration.z; 
                        input_in.acc = input_in.acc * G_m_s2 / acc_norm;
                        imu_deque.pop_front();
                        if (imu_deque.empty()) break;
                        imu_last = imu_next;
                        imu_next = *(imu_deque.front());
                        }
                        else
                        {
                            imu_deque.pop_front();
                            if (imu_deque.empty()) break;
                            imu_last = imu_next;
                            imu_next = *(imu_deque.front());
                        }
                    }
                    }
                }
            }
            // M3D rot_cur_lidar;
            // if (p_imu->UseLIInit)
            // {
            //     rot_cur_lidar = state.rot_end * Lidar_R_wrt_IMU;
            // }
            // else
            // {
            //     rot_cur_lidar = state.rot_end;
            // }
            // euler_cur = RotMtoEuler(rot_cur_lidar);
            // geoQuat = tf::createQuaternionMsgFromRollPitchYaw
            //                     (euler_cur(0), euler_cur(1), euler_cur(2));
            /******* Publish odometry downsample *******/
            if (!publish_odometry_without_downsample)
            {
                // publish_odometry(pubOdomAftMapped);
            }

            /*** add the feature points to map kdtree ***/
            t3 = omp_get_wtime();
            
            if(feats_down_size > 4)
            {
                map_incremental();
            }

            t5 = omp_get_wtime();
            /******* Publish points *******/
            if (path_en)                         publish_path(pubPath);
            if (scan_pub_en || pcd_save_en)      publish_frame_world(pubLaserCloudFullRes);
            if (scan_pub_en && scan_body_pub_en) publish_frame_body(pubLaserCloudFullRes_body);
            
            /*** Debug variables Logging ***/
            if (runtime_pos_log)
            {
                frame_num ++;
                aver_time_consu = aver_time_consu * (frame_num - 1) / frame_num + (t5 - t0) / frame_num;
                {aver_time_icp = aver_time_icp * (frame_num - 1)/frame_num + update_time/frame_num;}
                aver_time_match = aver_time_match * (frame_num - 1)/frame_num + (match_time)/frame_num;
                aver_time_solve = aver_time_solve * (frame_num - 1)/frame_num + solve_time/frame_num;
                aver_time_propag = aver_time_propag * (frame_num - 1)/frame_num + propag_time / frame_num;
                T1[time_log_counter] = Measures.lidar_beg_time;
                s_plot[time_log_counter] = t5 - t0;
                s_plot2[time_log_counter] = feats_undistort->points.size();
                s_plot3[time_log_counter] = aver_time_consu;
                time_log_counter ++;
                printf("[ mapping ]: time: IMU + Map + Input Downsample: %0.6f ave match: %0.6f ave solve: %0.6f  ave ICP: %0.6f  map incre: %0.6f ave total: %0.6f icp: %0.6f propogate: %0.6f \n",t1-t0,aver_time_match,aver_time_solve,t3-t1,t5-t3,aver_time_consu, aver_time_icp, aver_time_propag); 
                if (!publish_odometry_without_downsample)
                {
                    if (!use_imu_as_input)
                    {
                        if (!GNSS_ENABLE)
                        {
                            Eigen::Matrix3d R_enu_local_;
                            Eigen::Vector3d pos_r = kf_output.x_.rot.normalized() * p_gnss->Tex_imu_r + kf_output.x_.pos;
                            Eigen::Vector3d anc;
                            anc << offline_init_vec[0], offline_init_vec[1], offline_init_vec[2];
                            R_enu_local_ = ecef2rotation(anc) * Eigen::AngleAxisd(offline_init_vec[3], Eigen::Vector3d::UnitZ()); 
                            Eigen::Vector3d pos_enu = p_gnss->local2enu(R_enu_local_, anc, pos_r);
                            euler_cur = SO3ToEuler(kf_output.x_.rot);
                            fout_out << setw(20) << imu_next.header.stamp.toSec() - first_imu_time << " " << euler_cur.transpose() << " " << pos_enu.transpose() << " " << kf_output.x_.vel.transpose() \
                                        <<" "<<kf_output.x_.omg.transpose()<<" "<<kf_output.x_.acc.transpose()<<" "<<kf_output.x_.gravity.transpose()<<" "<<kf_output.x_.bg.transpose() <<" "<<kf_output.x_.ba.transpose()<<" "<< 0 <<endl;
                        }
                        // fout_out << setw(20) << Measures.lidar_beg_time - first_lidar_time << " " << euler_cur.transpose() << " " << kf_output.x_.pos.transpose() << " " << kf_output.x_.vel.transpose() \
                        // <<" "<<state_out.omg.transpose()<<" "<<state_out.acc.transpose()<<" "<<state_out.gravity.transpose()<<" "<<state_out.bg.transpose()<<" "<<state_out.ba.transpose()<<" "<<feats_undistort->points.size()<<endl;
                    }
                    else
                    {
                        if (!GNSS_ENABLE)
                        {
                            Eigen::Matrix3d R_enu_local_;
                            Eigen::Vector3d pos_r = kf_input.x_.rot.normalized() * p_gnss->Tex_imu_r + kf_input.x_.pos;
                            Eigen::Vector3d anc;
                            anc << offline_init_vec[0], offline_init_vec[1], offline_init_vec[2];
                            R_enu_local_ = ecef2rotation(anc) * Eigen::AngleAxisd(offline_init_vec[3], Eigen::Vector3d::UnitZ()); 
                            Eigen::Vector3d pos_enu = p_gnss->local2enu(R_enu_local_, anc, pos_r);
                            euler_cur = SO3ToEuler(kf_input.x_.rot);
                            fout_out << setw(20) << imu_next.header.stamp.toSec() - first_imu_time << " " << euler_cur.transpose() << " " << kf_input.x_.pos.transpose() << " " << kf_input.x_.vel.transpose() \
                            <<" "<<kf_input.x_.bg.transpose()<<" "<<kf_input.x_.ba.transpose()<<" "<<kf_input.x_.gravity.transpose()<<" " <<0<<endl;
                        }
                        // fout_out << setw(20) << Measures.lidar_beg_time - first_lidar_time << " " << euler_cur.transpose() << " " << state_in.pos.transpose() << " " << state_in.vel.transpose() \
                        // <<" "<<state_in.bg.transpose()<<" "<<state_in.ba.transpose()<<" "<<state_in.gravity.transpose()<<" "<<feats_undistort->points.size()<<endl;
                    }
                }
                dump_lio_state_to_log(fp);
            }
        }
        status = ros::ok();
        rate.sleep();
    }
    //--------------------------save map-----------------------------------
    /* 1. make sure you have enough memories
    /* 2. noted that pcd save will influence the real-time performences **/
    if (pcl_wait_save->size() > 0 && pcd_save_en)
    {
        string file_name = string("scans.pcd");
        string all_points_dir(string(string(ROOT_DIR) + "PCD/") + file_name);
        pcl::PCDWriter pcd_writer;
        pcd_writer.writeBinary(all_points_dir, *pcl_wait_save);
    }
    fout_out.close();
    fout_imu_pbp.close();
    for (int i = 0; i < p_gnss->pvt_time.size(); i++)
    {
        fout_rtk << setw(20) << p_gnss->pvt_time[i] - p_gnss->pvt_time[0] << " " << p_gnss->pvt_holder[i].transpose() << endl;
    }
    fout_rtk.close();
    return 0;
}
