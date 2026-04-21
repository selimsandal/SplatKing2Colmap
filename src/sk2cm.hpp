#pragma once

#include <cstdint>
#include <filesystem>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace sk2cm {

namespace fs = std::filesystem;

struct Vector3d
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vector3d operator+(const Vector3d& lhs, const Vector3d& rhs);
Vector3d operator-(const Vector3d& lhs, const Vector3d& rhs);
Vector3d operator*(const Vector3d& value, double scalar);
Vector3d operator*(double scalar, const Vector3d& value);
Vector3d operator/(const Vector3d& value, double scalar);
Vector3d& operator+=(Vector3d& lhs, const Vector3d& rhs);
double Dot(const Vector3d& lhs, const Vector3d& rhs);
double Norm(const Vector3d& value);
double Distance(const Vector3d& lhs, const Vector3d& rhs);
double DistanceSquared(const Vector3d& lhs, const Vector3d& rhs);
Vector3d Normalize(const Vector3d& value);

struct Matrix3d
{
    double m11 = 1.0;
    double m12 = 0.0;
    double m13 = 0.0;
    double m21 = 0.0;
    double m22 = 1.0;
    double m23 = 0.0;
    double m31 = 0.0;
    double m32 = 0.0;
    double m33 = 1.0;
};

Matrix3d IdentityMatrix3d();
Matrix3d Transpose(const Matrix3d& value);
Matrix3d Multiply(const Matrix3d& lhs, const Matrix3d& rhs);
Vector3d Multiply(const Matrix3d& lhs, const Vector3d& rhs);

struct Quaterniond
{
    double w = 1.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Quaterniond Normalize(const Quaterniond& value);
Quaterniond Negate(const Quaterniond& value);
Matrix3d QuaternionToRotationMatrix(const Quaterniond& value);
Quaterniond RotationMatrixToHamiltonQuaternion(const Matrix3d& rotation);

struct ColmapPose
{
    Matrix3d world_to_camera;
    Quaterniond rotation;
    Vector3d translation;
};

struct DepthMapInfo
{
    fs::path path;
    int width = 0;
    int height = 0;
};

enum class ImageRotation
{
    None = 0,
    Rotate90Clockwise = 1,
    Rotate90CounterClockwise = 2
};

struct LoadedCapture
{
    int index = 0;
    std::string image_name;
    fs::path image_path;
    fs::path metadata_path;
    std::optional<fs::path> depth_path;
    int width = 0;
    int height = 0;
    double timestamp = 0.0;
    double fx = 0.0;
    double fy = 0.0;
    double cx = 0.0;
    double cy = 0.0;
    double metadata_fx = 0.0;
    double metadata_fy = 0.0;
    double metadata_cx = 0.0;
    double metadata_cy = 0.0;
    int metadata_width = 0;
    int metadata_height = 0;
    std::string tracking_state;
    ImageRotation image_rotation = ImageRotation::None;
    Vector3d camera_center_world;
    ColmapPose metadata_colmap_pose;
    ColmapPose colmap_pose;
    std::optional<DepthMapInfo> depth_info;
};

struct LoadedDataset
{
    fs::path input_directory;
    std::vector<LoadedCapture> captures;
    std::vector<Vector3d> world_points;
    fs::path point_cloud_path;
};

struct ProjectionSample
{
    double pixel_x = 0.0;
    double pixel_y = 0.0;
    double depth_meters = 0.0;
};

struct ValidationSummary
{
    std::string label;
    int images_checked = 0;
    int points_checked = 0;
    int positive_depth_count = 0;
    int in_frame_count = 0;
    int depth_compared_count = 0;
    double mean_depth_abs_error = 0.0;
    double max_depth_abs_error = 0.0;
};

struct ReprojectionSummary
{
    int point_count = 0;
    int observation_count = 0;
    double mean_pixel_error = 0.0;
    double max_pixel_error = 0.0;
};

struct CameraPoseSummary
{
    int camera_count = 0;
    double mean_center_error_meters = 0.0;
    double max_center_error_meters = 0.0;
};

struct AlignmentSummary
{
    int matched_image_count = 0;
    int snapped_image_count = 0;
    bool shared_camera_poses_snapped_to_reference = false;
    double scale = 1.0;
    double mean_center_residual_meters = 0.0;
    double max_center_residual_meters = 0.0;
    std::string reference_model_path;
};

struct LidarObservationCandidate
{
    int point_index = 0;
    double pixel_x = 0.0;
    double pixel_y = 0.0;
    double depth_meters = 0.0;
    double depth_abs_error = 0.0;
    int depth_cell_index = 0;
};

struct PointObservationCandidate
{
    int image_id = 0;
    int capture_index = 0;
    double pixel_x = 0.0;
    double pixel_y = 0.0;
    double depth_abs_error = 0.0;
    Vector3d camera_center_world;
};

struct ImageObservation
{
    double x = 0.0;
    double y = 0.0;
    std::int64_t point3d_id = 0;
};

struct PointTrackElement
{
    int image_id = 0;
    int point2d_index = 0;
};

struct SparseImagePoints
{
    int image_id = 0;
    std::vector<ImageObservation> observations;
};

struct SparsePoint3D
{
    std::int64_t point3d_id = 0;
    Vector3d position;
    std::uint8_t r = 255;
    std::uint8_t g = 255;
    std::uint8_t b = 255;
    double error = 0.0;
    std::vector<PointTrackElement> track;
};

struct SparseReconstruction
{
    std::vector<SparseImagePoints> image_points;
    std::vector<SparsePoint3D> points3d;
};

struct ReconstructionBuildResult
{
    SparseReconstruction sparse_reconstruction;
    std::vector<Vector3d> dense_aligned_points;
};

struct AlignmentResult
{
    LoadedDataset dataset;
    ReconstructionBuildResult reconstruction;
    AlignmentSummary summary;
};

struct LidarRefinementOptions
{
    int depth_bucket_size = 4;
    int max_observations_per_track = 6;
    int min_capture_gap = 4;
    double min_baseline_meters = 0.12;
    double min_angular_separation_degrees = 2.0;
    double base_depth_tolerance_meters = 0.50;
    double relative_depth_tolerance = 0.20;
};

struct ColmapModelCounts
{
    std::uint64_t cameras = 0;
    std::uint64_t images = 0;
    std::uint64_t points3d = 0;
    std::uint64_t image_observations = 0;
    std::uint64_t point_track_elements = 0;
};

struct ColmapModelValidationSummary
{
    ColmapModelCounts text_counts;
    ColmapModelCounts binary_counts;

    bool IsMatch() const;
};

struct CliOptions
{
    fs::path input_directory;
    fs::path output_directory;
    std::optional<fs::path> align_to_colmap_model_path;
    bool snap_shared_camera_poses_to_reference = false;
    bool copy_images = true;
    bool validate_exported_model = false;
    int validation_image_count = 5;
    int validation_point_stride = 256;
    LidarRefinementOptions refinement;
};

struct ConversionResult
{
    fs::path output_directory;
    int capture_count = 0;
    int dense_point_count = 0;
    int sparse_point_count = 0;
    int observation_count = 0;
    double average_track_length = 0.0;
    std::optional<AlignmentSummary> alignment;
    ValidationSummary validation;
    std::optional<ColmapModelValidationSummary> model_validation;
};

class HelpRequested final : public std::runtime_error
{
public:
    HelpRequested();
};

CliOptions CreateDefaultOptions(
    const fs::path& input_directory,
    const fs::path& output_directory,
    bool validate_exported_model = true,
    bool copy_images = true);

CliOptions ParseCliOptions(const std::vector<std::string>& args);
std::string GetUsageText();

ConversionResult RunConversion(
    const CliOptions& options,
    const std::function<void(const std::string&)>& log = {});

}  // namespace sk2cm
