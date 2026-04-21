#include "sk2cm.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <ios>
#include <istream>
#include <iterator>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <ostream>
#include <set>
#include <sstream>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include <nlohmann/json.hpp>

namespace sk2cm {

Vector3d operator+(const Vector3d& lhs, const Vector3d& rhs)
{
    return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Vector3d operator-(const Vector3d& lhs, const Vector3d& rhs)
{
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Vector3d operator*(const Vector3d& value, double scalar)
{
    return {value.x * scalar, value.y * scalar, value.z * scalar};
}

Vector3d operator*(double scalar, const Vector3d& value)
{
    return value * scalar;
}

Vector3d operator/(const Vector3d& value, double scalar)
{
    return {value.x / scalar, value.y / scalar, value.z / scalar};
}

Vector3d& operator+=(Vector3d& lhs, const Vector3d& rhs)
{
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    return lhs;
}

double Dot(const Vector3d& lhs, const Vector3d& rhs)
{
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

double Norm(const Vector3d& value)
{
    return std::sqrt(Dot(value, value));
}

double Distance(const Vector3d& lhs, const Vector3d& rhs)
{
    return Norm(lhs - rhs);
}

double DistanceSquared(const Vector3d& lhs, const Vector3d& rhs)
{
    const auto delta = lhs - rhs;
    return Dot(delta, delta);
}

Vector3d Normalize(const Vector3d& value)
{
    const auto norm = Norm(value);
    if (norm <= 1e-12)
    {
        return {};
    }

    return value / norm;
}

Matrix3d IdentityMatrix3d()
{
    return {};
}

Matrix3d Transpose(const Matrix3d& value)
{
    return {
        value.m11, value.m21, value.m31,
        value.m12, value.m22, value.m32,
        value.m13, value.m23, value.m33};
}

Matrix3d Multiply(const Matrix3d& lhs, const Matrix3d& rhs)
{
    return {
        (lhs.m11 * rhs.m11) + (lhs.m12 * rhs.m21) + (lhs.m13 * rhs.m31),
        (lhs.m11 * rhs.m12) + (lhs.m12 * rhs.m22) + (lhs.m13 * rhs.m32),
        (lhs.m11 * rhs.m13) + (lhs.m12 * rhs.m23) + (lhs.m13 * rhs.m33),
        (lhs.m21 * rhs.m11) + (lhs.m22 * rhs.m21) + (lhs.m23 * rhs.m31),
        (lhs.m21 * rhs.m12) + (lhs.m22 * rhs.m22) + (lhs.m23 * rhs.m32),
        (lhs.m21 * rhs.m13) + (lhs.m22 * rhs.m23) + (lhs.m23 * rhs.m33),
        (lhs.m31 * rhs.m11) + (lhs.m32 * rhs.m21) + (lhs.m33 * rhs.m31),
        (lhs.m31 * rhs.m12) + (lhs.m32 * rhs.m22) + (lhs.m33 * rhs.m32),
        (lhs.m31 * rhs.m13) + (lhs.m32 * rhs.m23) + (lhs.m33 * rhs.m33)};
}

Vector3d Multiply(const Matrix3d& lhs, const Vector3d& rhs)
{
    return {
        (lhs.m11 * rhs.x) + (lhs.m12 * rhs.y) + (lhs.m13 * rhs.z),
        (lhs.m21 * rhs.x) + (lhs.m22 * rhs.y) + (lhs.m23 * rhs.z),
        (lhs.m31 * rhs.x) + (lhs.m32 * rhs.y) + (lhs.m33 * rhs.z)};
}

Quaterniond Normalize(const Quaterniond& value)
{
    const auto norm = std::sqrt(
        (value.w * value.w) +
        (value.x * value.x) +
        (value.y * value.y) +
        (value.z * value.z));
    if (norm <= 1e-12)
    {
        return {};
    }

    return {
        value.w / norm,
        value.x / norm,
        value.y / norm,
        value.z / norm};
}

Quaterniond Negate(const Quaterniond& value)
{
    return {-value.w, -value.x, -value.y, -value.z};
}

Matrix3d QuaternionToRotationMatrix(const Quaterniond& input)
{
    const auto q = Normalize(input);
    return {
        1.0 - (2.0 * q.y * q.y) - (2.0 * q.z * q.z),
        (2.0 * q.x * q.y) - (2.0 * q.z * q.w),
        (2.0 * q.x * q.z) + (2.0 * q.y * q.w),
        (2.0 * q.x * q.y) + (2.0 * q.z * q.w),
        1.0 - (2.0 * q.x * q.x) - (2.0 * q.z * q.z),
        (2.0 * q.y * q.z) - (2.0 * q.x * q.w),
        (2.0 * q.x * q.z) - (2.0 * q.y * q.w),
        (2.0 * q.y * q.z) + (2.0 * q.x * q.w),
        1.0 - (2.0 * q.x * q.x) - (2.0 * q.y * q.y)};
}

Quaterniond RotationMatrixToHamiltonQuaternion(const Matrix3d& rotation)
{
    const auto rxx = rotation.m11;
    const auto rxy = rotation.m12;
    const auto rxz = rotation.m13;
    const auto ryx = rotation.m21;
    const auto ryy = rotation.m22;
    const auto ryz = rotation.m23;
    const auto rzx = rotation.m31;
    const auto rzy = rotation.m32;
    const auto rzz = rotation.m33;

    double qw = 0.0;
    double qx = 0.0;
    double qy = 0.0;
    double qz = 0.0;

    const auto trace = rxx + ryy + rzz;
    if (trace > 0.0)
    {
        const auto s = std::sqrt(trace + 1.0) * 2.0;
        qw = 0.25 * s;
        qx = (rzy - ryz) / s;
        qy = (rxz - rzx) / s;
        qz = (ryx - rxy) / s;
    }
    else if (rxx > ryy && rxx > rzz)
    {
        const auto s = std::sqrt(1.0 + rxx - ryy - rzz) * 2.0;
        qw = (rzy - ryz) / s;
        qx = 0.25 * s;
        qy = (rxy + ryx) / s;
        qz = (rxz + rzx) / s;
    }
    else if (ryy > rzz)
    {
        const auto s = std::sqrt(1.0 + ryy - rxx - rzz) * 2.0;
        qw = (rxz - rzx) / s;
        qx = (rxy + ryx) / s;
        qy = 0.25 * s;
        qz = (ryz + rzy) / s;
    }
    else
    {
        const auto s = std::sqrt(1.0 + rzz - rxx - ryy) * 2.0;
        qw = (ryx - rxy) / s;
        qx = (rxz + rzx) / s;
        qy = (ryz + rzy) / s;
        qz = 0.25 * s;
    }

    auto quaternion = Normalize(Quaterniond{qw, qx, qy, qz});
    if (quaternion.w < 0.0)
    {
        quaternion = Negate(quaternion);
    }

    return quaternion;
}

namespace {

using json = nlohmann::json;

struct ParsedArkitTransform
{
    Matrix3d rotation;
    Vector3d translation;
};

struct DepthPixel
{
    int x = 0;
    int y = 0;
};

struct ReferenceImagePose
{
    std::string image_name;
    Vector3d camera_center_world;
    ColmapPose pose;
};

struct SimilarityTransform
{
    Matrix3d rotation;
    Vector3d translation;
    double scale = 1.0;
};

struct VoxelKey
{
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const VoxelKey& other) const = default;
};

struct VoxelKeyHash
{
    std::size_t operator()(const VoxelKey& key) const noexcept
    {
        auto seed = static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.x));
        seed ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.y)) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        seed ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.z)) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
        return static_cast<std::size_t>(seed);
    }
};

class TargetPointIndex
{
public:
    static TargetPointIndex Build(const std::vector<Vector3d>& points, double voxel_size_meters)
    {
        std::unordered_map<VoxelKey, std::vector<Vector3d>, VoxelKeyHash> points_by_voxel;
        points_by_voxel.reserve(points.size());
        for (const auto& point : points)
        {
            auto key = FromPoint(point, voxel_size_meters);
            points_by_voxel[key].push_back(point);
        }

        return TargetPointIndex(voxel_size_meters, std::move(points_by_voxel));
    }

    bool TryFindNearest(const Vector3d& point, double max_distance_meters, Vector3d& nearest) const
    {
        const auto center_key = FromPoint(point, voxel_size_meters_);
        const auto radius_in_cells = std::max(1, static_cast<int>(std::ceil(max_distance_meters / voxel_size_meters_)));
        auto best_squared_distance = max_distance_meters * max_distance_meters;
        auto found = false;

        for (int dz = -radius_in_cells; dz <= radius_in_cells; ++dz)
        {
            for (int dy = -radius_in_cells; dy <= radius_in_cells; ++dy)
            {
                for (int dx = -radius_in_cells; dx <= radius_in_cells; ++dx)
                {
                    const VoxelKey key{center_key.x + dx, center_key.y + dy, center_key.z + dz};
                    const auto it = points_by_voxel_.find(key);
                    if (it == points_by_voxel_.end())
                    {
                        continue;
                    }

                    for (const auto& candidate : it->second)
                    {
                        const auto squared_distance = DistanceSquared(point, candidate);
                        if (squared_distance < best_squared_distance)
                        {
                            best_squared_distance = squared_distance;
                            nearest = candidate;
                            found = true;
                        }
                    }
                }
            }
        }

        return found;
    }

private:
    TargetPointIndex(
        double voxel_size_meters,
        std::unordered_map<VoxelKey, std::vector<Vector3d>, VoxelKeyHash> points_by_voxel)
        : voxel_size_meters_(voxel_size_meters),
          points_by_voxel_(std::move(points_by_voxel))
    {
    }

    static VoxelKey FromPoint(const Vector3d& point, double voxel_size_meters)
    {
        return {
            static_cast<int>(std::floor(point.x / voxel_size_meters)),
            static_cast<int>(std::floor(point.y / voxel_size_meters)),
            static_cast<int>(std::floor(point.z / voxel_size_meters))};
    }

    double voxel_size_meters_ = 1.0;
    std::unordered_map<VoxelKey, std::vector<Vector3d>, VoxelKeyHash> points_by_voxel_;
};

std::string BoolText(const bool value)
{
    return value ? "True" : "False";
}

std::string ToRoundTrip(const double value)
{
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<double>::max_digits10) << value;
    return stream.str();
}

std::string ToFixed(const double value, const int decimals)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(decimals) << value;
    return stream.str();
}

std::string ReadTextFile(const fs::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + path.string() + "'.");
    }

    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

std::vector<std::string> SplitWhitespace(const std::string& line)
{
    std::vector<std::string> parts;
    std::istringstream stream(line);
    std::string part;
    while (stream >> part)
    {
        parts.push_back(part);
    }

    return parts;
}

bool IEquals(std::string_view lhs, std::string_view rhs)
{
    if (lhs.size() != rhs.size())
    {
        return false;
    }

    for (std::size_t index = 0; index < lhs.size(); ++index)
    {
        const auto l = static_cast<unsigned char>(lhs[index]);
        const auto r = static_cast<unsigned char>(rhs[index]);
        if (std::tolower(l) != std::tolower(r))
        {
            return false;
        }
    }

    return true;
}

int Clamp(const int value, const int minimum, const int maximum)
{
    return std::clamp(value, minimum, maximum);
}

double RoundToEven(const double value)
{
    const auto floor_value = std::floor(value);
    const auto fraction = value - floor_value;
    if (fraction < 0.5)
    {
        return floor_value;
    }

    if (fraction > 0.5)
    {
        return floor_value + 1.0;
    }

    const auto floor_as_int = static_cast<long long>(floor_value);
    return (floor_as_int % 2LL) == 0LL ? floor_value : (floor_value + 1.0);
}

template <typename T>
void WriteLittleEndian(std::ostream& stream, const T value)
{
    static_assert(std::is_trivially_copyable_v<T>);
    auto bytes = std::bit_cast<std::array<std::byte, sizeof(T)>>(value);
    if constexpr (std::endian::native == std::endian::big)
    {
        std::reverse(bytes.begin(), bytes.end());
    }

    stream.write(reinterpret_cast<const char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!stream)
    {
        throw std::runtime_error("Failed to write binary stream.");
    }
}

template <typename T>
T ReadLittleEndian(std::istream& stream)
{
    static_assert(std::is_trivially_copyable_v<T>);
    std::array<std::byte, sizeof(T)> bytes{};
    stream.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!stream)
    {
        throw std::runtime_error("Failed to read binary stream.");
    }

    if constexpr (std::endian::native == std::endian::big)
    {
        std::reverse(bytes.begin(), bytes.end());
    }

    return std::bit_cast<T>(bytes);
}

double ParseDoubleInvariant(const std::string& value, const std::string& option_name)
{
    std::size_t consumed = 0;
    const auto result = std::stod(value, &consumed);
    if (consumed != value.size())
    {
        throw std::invalid_argument("Option '" + option_name + "' expects a number, got '" + value + "'.");
    }

    return result;
}

int ParsePositiveInt(const std::string& value, const std::string& option_name)
{
    std::size_t consumed = 0;
    const auto result = std::stoi(value, &consumed, 10);
    if (consumed != value.size() || result <= 0)
    {
        throw std::invalid_argument("Option '" + option_name + "' expects a positive integer, got '" + value + "'.");
    }

    return result;
}

int ParseNonNegativeInt(const std::string& value, const std::string& option_name)
{
    std::size_t consumed = 0;
    const auto result = std::stoi(value, &consumed, 10);
    if (consumed != value.size() || result < 0)
    {
        throw std::invalid_argument("Option '" + option_name + "' expects a non-negative integer, got '" + value + "'.");
    }

    return result;
}

double ParsePositiveDouble(const std::string& value, const std::string& option_name)
{
    const auto result = ParseDoubleInvariant(value, option_name);
    if (!(result > 0.0))
    {
        throw std::invalid_argument("Option '" + option_name + "' expects a positive number, got '" + value + "'.");
    }

    return result;
}

std::string RequireValue(const std::vector<std::string>& args, std::size_t& index, const std::string& option_name)
{
    if (index + 1 >= args.size())
    {
        throw std::invalid_argument("Option '" + option_name + "' requires a value.");
    }

    ++index;
    return args[index];
}

std::pair<int, int> ReadJpegDimensions(const fs::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open JPEG '" + path.string() + "'.");
    }

    const auto first = ReadLittleEndian<std::uint8_t>(stream);
    const auto second = ReadLittleEndian<std::uint8_t>(stream);
    if (first != 0xFFU || second != 0xD8U)
    {
        throw std::runtime_error("File '" + path.string() + "' is not a valid JPEG.");
    }

    auto read_big_endian_u16 = [&stream]() -> std::uint16_t {
        const auto high = ReadLittleEndian<std::uint8_t>(stream);
        const auto low = ReadLittleEndian<std::uint8_t>(stream);
        return static_cast<std::uint16_t>((static_cast<std::uint16_t>(high) << 8U) | static_cast<std::uint16_t>(low));
    };

    while (stream.good())
    {
        std::uint8_t marker_prefix = 0;
        do
        {
            stream.read(reinterpret_cast<char*>(&marker_prefix), 1);
            if (!stream)
            {
                break;
            }
        } while (marker_prefix != 0xFFU);

        if (!stream)
        {
            break;
        }

        std::uint8_t marker = 0;
        do
        {
            stream.read(reinterpret_cast<char*>(&marker), 1);
            if (!stream)
            {
                break;
            }
        } while (marker == 0xFFU);

        if (!stream)
        {
            break;
        }

        if (marker == 0xD8U || marker == 0xD9U)
        {
            continue;
        }

        const auto segment_length = read_big_endian_u16();
        if (segment_length < 2U)
        {
            throw std::runtime_error("Invalid JPEG segment length in '" + path.string() + "'.");
        }

        switch (marker)
        {
            case 0xC0U:
            case 0xC1U:
            case 0xC2U:
            case 0xC3U:
            case 0xC5U:
            case 0xC6U:
            case 0xC7U:
            case 0xC9U:
            case 0xCAU:
            case 0xCBU:
            case 0xCDU:
            case 0xCEU:
            case 0xCFU:
            {
                stream.ignore(1);
                const auto height = static_cast<int>(read_big_endian_u16());
                const auto width = static_cast<int>(read_big_endian_u16());
                return {width, height};
            }
            default:
                stream.seekg(static_cast<std::streamoff>(segment_length) - 2, std::ios::cur);
                break;
        }
    }

    throw std::runtime_error("Could not read JPEG dimensions from '" + path.string() + "'.");
}

std::vector<float> LoadDepthMap(const DepthMapInfo& depth_info)
{
    std::ifstream stream(depth_info.path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open depth file '" + depth_info.path.string() + "'.");
    }

    const auto count = static_cast<std::size_t>(depth_info.width) * static_cast<std::size_t>(depth_info.height);
    const auto expected_bytes = count * sizeof(float);
    stream.seekg(0, std::ios::end);
    const auto actual_bytes = static_cast<std::size_t>(stream.tellg());
    stream.seekg(0, std::ios::beg);
    if (actual_bytes != expected_bytes)
    {
        throw std::runtime_error(
            "Depth file '" + depth_info.path.string() + "' has " +
            std::to_string(actual_bytes) +
            " bytes, expected " +
            std::to_string(expected_bytes) +
            ".");
    }

    std::vector<float> values(count);
    for (auto& value : values)
    {
        value = ReadLittleEndian<float>(stream);
    }

    return values;
}

std::vector<Vector3d> LoadFloat32PointCloud(const fs::path& point_cloud_path)
{
    std::ifstream stream(point_cloud_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Missing LiDAR point cloud binary '" + point_cloud_path.string() + "'.");
    }

    stream.seekg(0, std::ios::end);
    const auto total_bytes = static_cast<std::size_t>(stream.tellg());
    stream.seekg(0, std::ios::beg);
    if ((total_bytes % 12U) != 0U)
    {
        throw std::runtime_error("Point cloud binary length is not divisible by 12 bytes.");
    }

    std::vector<Vector3d> points(total_bytes / 12U);
    for (auto& point : points)
    {
        point.x = static_cast<double>(ReadLittleEndian<float>(stream));
        point.y = static_cast<double>(ReadLittleEndian<float>(stream));
        point.z = static_cast<double>(ReadLittleEndian<float>(stream));
    }

    return points;
}

ParsedArkitTransform ParseArkitCameraToWorld(const std::vector<double>& values)
{
    if (values.size() != 16U)
    {
        throw std::runtime_error("Expected 16 transform values, got " + std::to_string(values.size()) + ".");
    }

    return {
        {
            values[0], values[4], values[8],
            values[1], values[5], values[9],
            values[2], values[6], values[10]},
        {values[12], values[13], values[14]}};
}

ColmapPose ConvertArkitCameraToColmap(const ParsedArkitTransform& arkit_camera_to_world)
{
    const auto& row0 = Vector3d{
        arkit_camera_to_world.rotation.m11,
        arkit_camera_to_world.rotation.m12,
        arkit_camera_to_world.rotation.m13};
    const auto& row1 = Vector3d{
        arkit_camera_to_world.rotation.m21,
        arkit_camera_to_world.rotation.m22,
        arkit_camera_to_world.rotation.m23};
    const auto& row2 = Vector3d{
        arkit_camera_to_world.rotation.m31,
        arkit_camera_to_world.rotation.m32,
        arkit_camera_to_world.rotation.m33};

    const auto world_to_camera_rotation = Matrix3d{
        row0.x, row1.x, row2.x,
        -row0.y, -row1.y, -row2.y,
        -row0.z, -row1.z, -row2.z};

    const auto translation_world = arkit_camera_to_world.translation;
    const auto translation = Vector3d{
        -Dot({world_to_camera_rotation.m11, world_to_camera_rotation.m12, world_to_camera_rotation.m13}, translation_world),
        -Dot({world_to_camera_rotation.m21, world_to_camera_rotation.m22, world_to_camera_rotation.m23}, translation_world),
        -Dot({world_to_camera_rotation.m31, world_to_camera_rotation.m32, world_to_camera_rotation.m33}, translation_world)};

    return {
        world_to_camera_rotation,
        RotationMatrixToHamiltonQuaternion(world_to_camera_rotation),
        translation};
}

ColmapPose RotateImage90Clockwise(const ColmapPose& pose)
{
    const auto image_rotation = Matrix3d{
        0.0, -1.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 1.0};
    const auto rotated_world_to_camera = Multiply(image_rotation, pose.world_to_camera);
    const auto rotated_translation = Multiply(image_rotation, pose.translation);
    return {
        rotated_world_to_camera,
        RotationMatrixToHamiltonQuaternion(rotated_world_to_camera),
        rotated_translation};
}

ColmapPose RotateImage90CounterClockwise(const ColmapPose& pose)
{
    const auto image_rotation = Matrix3d{
        0.0, 1.0, 0.0,
        -1.0, 0.0, 0.0,
        0.0, 0.0, 1.0};
    const auto rotated_world_to_camera = Multiply(image_rotation, pose.world_to_camera);
    const auto rotated_translation = Multiply(image_rotation, pose.translation);
    return {
        rotated_world_to_camera,
        RotationMatrixToHamiltonQuaternion(rotated_world_to_camera),
        rotated_translation};
}

std::optional<ProjectionSample> ProjectWorldPoint(const LoadedCapture& capture, const Vector3d& world_point)
{
    const auto camera_point = Multiply(capture.colmap_pose.world_to_camera, world_point) + capture.colmap_pose.translation;
    if (camera_point.z <= 1e-6)
    {
        return std::nullopt;
    }

    return ProjectionSample{
        capture.fx * (camera_point.x / camera_point.z) + capture.cx,
        capture.fy * (camera_point.y / camera_point.z) + capture.cy,
        camera_point.z};
}

std::optional<std::pair<double, double>> MapImagePixelToMetadataPixel(
    const LoadedCapture& capture,
    const double pixel_x,
    const double pixel_y)
{
    switch (capture.image_rotation)
    {
        case ImageRotation::None:
            return std::make_pair(pixel_x, pixel_y);
        case ImageRotation::Rotate90Clockwise:
            return std::make_pair(pixel_y, static_cast<double>(capture.metadata_height) - pixel_x);
        case ImageRotation::Rotate90CounterClockwise:
            return std::make_pair(static_cast<double>(capture.metadata_width) - pixel_y, pixel_x);
        default:
            throw std::runtime_error("Unsupported image rotation.");
    }
}

std::optional<DepthPixel> MapImagePixelToDepthPixel(
    const LoadedCapture& capture,
    const double pixel_x,
    const double pixel_y)
{
    if (!capture.depth_info.has_value())
    {
        return std::nullopt;
    }

    const auto metadata_pixel = MapImagePixelToMetadataPixel(capture, pixel_x, pixel_y);
    if (!metadata_pixel.has_value())
    {
        return std::nullopt;
    }

    const auto& depth_info = *capture.depth_info;
    const auto depth_x = Clamp(
        static_cast<int>(RoundToEven(
            metadata_pixel->first * static_cast<double>(depth_info.width) / static_cast<double>(capture.metadata_width))),
        0,
        depth_info.width - 1);
    const auto depth_y = Clamp(
        static_cast<int>(RoundToEven(
            metadata_pixel->second * static_cast<double>(depth_info.height) / static_cast<double>(capture.metadata_height))),
        0,
        depth_info.height - 1);

    return DepthPixel{depth_x, depth_y};
}

LoadedCapture LoadCapture(const fs::path& input_directory, const json& capture_json, const int index)
{
    const auto image_name = capture_json.at("filename").get<std::string>();
    const auto image_path = input_directory / image_name;
    if (!fs::exists(image_path))
    {
        throw std::runtime_error("Missing image '" + image_name + "'.");
    }

    const auto [actual_image_width, actual_image_height] = ReadJpegDimensions(image_path);

    const auto metadata_file = capture_json.at("metadataFile").get<std::string>();
    const auto metadata_path = input_directory / metadata_file;
    if (!fs::exists(metadata_path))
    {
        throw std::runtime_error("Missing metadata file '" + metadata_file + "'.");
    }

    const auto& metadata_json = capture_json.at("metadata");
    const auto& camera_json = metadata_json.at("camera");
    const auto transform_values = camera_json.at("transform").get<std::vector<double>>();
    const auto intrinsics = camera_json.at("intrinsics").get<std::vector<double>>();
    if (transform_values.size() != 16U)
    {
        throw std::runtime_error("Capture '" + image_name + "' has an invalid transform size.");
    }

    if (intrinsics.size() != 9U)
    {
        throw std::runtime_error("Capture '" + image_name + "' has an invalid intrinsics size.");
    }

    std::optional<fs::path> depth_path;
    std::optional<DepthMapInfo> depth_info;
    if (capture_json.contains("auxiliaryOutputs") && capture_json.at("auxiliaryOutputs").is_array())
    {
        for (const auto& output_json : capture_json.at("auxiliaryOutputs"))
        {
            const auto type = output_json.value("type", std::string{});
            if (!IEquals(type, "depth"))
            {
                continue;
            }

            const auto filename = output_json.at("filename").get<std::string>();
            const auto candidate_path = input_directory / filename;
            if (!fs::exists(candidate_path))
            {
                throw std::runtime_error("Missing depth file '" + filename + "'.");
            }

            depth_path = candidate_path;
            depth_info = DepthMapInfo{
                candidate_path,
                output_json.at("width").get<int>(),
                output_json.at("height").get<int>()};
            break;
        }
    }

    const auto arkit_camera_to_world = ParseArkitCameraToWorld(transform_values);
    auto metadata_colmap_pose = ConvertArkitCameraToColmap(arkit_camera_to_world);
    auto colmap_pose = metadata_colmap_pose;
    auto width = camera_json.at("imageResolution").at("width").get<int>();
    auto height = camera_json.at("imageResolution").at("height").get<int>();
    const auto metadata_fx = intrinsics[0];
    const auto metadata_fy = intrinsics[4];
    const auto metadata_cx = intrinsics[6];
    const auto metadata_cy = intrinsics[7];
    auto fx = metadata_fx;
    auto fy = metadata_fy;
    auto cx = metadata_cx;
    auto cy = metadata_cy;
    auto image_rotation = ImageRotation::None;

    if (actual_image_width == height && actual_image_height == width)
    {
        colmap_pose = RotateImage90Clockwise(colmap_pose);
        const auto rotated_width = actual_image_width;
        const auto rotated_height = actual_image_height;
        const auto rotated_fx = fy;
        const auto rotated_fy = fx;
        const auto rotated_cx = static_cast<double>(height) - cy;
        const auto rotated_cy = cx;
        width = rotated_width;
        height = rotated_height;
        fx = rotated_fx;
        fy = rotated_fy;
        cx = rotated_cx;
        cy = rotated_cy;
        image_rotation = ImageRotation::Rotate90Clockwise;
    }
    else
    {
        width = actual_image_width;
        height = actual_image_height;
    }

    return {
        index,
        image_name,
        image_path,
        metadata_path,
        depth_path,
        width,
        height,
        capture_json.value("timestamp", 0.0),
        fx,
        fy,
        cx,
        cy,
        metadata_fx,
        metadata_fy,
        metadata_cx,
        metadata_cy,
        camera_json.at("imageResolution").at("width").get<int>(),
        camera_json.at("imageResolution").at("height").get<int>(),
        camera_json.value("trackingState", std::string{}),
        image_rotation,
        arkit_camera_to_world.translation,
        metadata_colmap_pose,
        colmap_pose,
        depth_info};
}

LoadedDataset LoadDataset(const fs::path& input_directory)
{
    const auto photo_series_path = input_directory / "photo_series.json";
    if (!fs::exists(photo_series_path))
    {
        throw std::runtime_error("Missing photo_series.json");
    }

    const auto photo_series = json::parse(ReadTextFile(photo_series_path));
    if (!photo_series.contains("captures") || !photo_series.at("captures").is_array())
    {
        throw std::runtime_error("Failed to deserialize photo_series.json.");
    }

    const auto& captures_json = photo_series.at("captures");
    if (captures_json.empty())
    {
        throw std::runtime_error("No captures were found in photo_series.json.");
    }

    std::string point_cloud_file_name;
    for (const auto& capture_json : captures_json)
    {
        const auto name = capture_json
            .value("extraMetadata", json::object())
            .value("sessionPointCloudFilename", std::string{});
        if (!name.empty())
        {
            point_cloud_file_name = name;
            break;
        }
    }

    if (point_cloud_file_name.empty())
    {
        throw std::runtime_error("Could not determine the LiDAR point cloud file name.");
    }

    const auto point_cloud_path = input_directory / point_cloud_file_name;
    std::vector<LoadedCapture> captures;
    captures.reserve(captures_json.size());
    for (std::size_t index = 0; index < captures_json.size(); ++index)
    {
        captures.push_back(LoadCapture(input_directory, captures_json[index], static_cast<int>(index + 1U)));
    }

    auto world_points = LoadFloat32PointCloud(point_cloud_path);
    const auto expected_point_count = captures_json.back()
        .value("extraMetadata", json::object())
        .value("sessionPointCloudPointCount", 0);
    if (expected_point_count > 0 && expected_point_count != static_cast<int>(world_points.size()))
    {
        throw std::runtime_error(
            "Point cloud count mismatch. JSON advertises " +
            std::to_string(expected_point_count) +
            " points but the binary contains " +
            std::to_string(world_points.size()) +
            ".");
    }

    return {
        fs::absolute(input_directory),
        std::move(captures),
        std::move(world_points),
        point_cloud_path};
}

bool HasGeometricSpread(
    const Vector3d& world_point,
    const PointObservationCandidate& first,
    const PointObservationCandidate& second,
    const LidarRefinementOptions& options)
{
    const auto baseline = Distance(first.camera_center_world, second.camera_center_world);
    if (baseline < options.min_baseline_meters)
    {
        return false;
    }

    const auto first_view = Normalize(world_point - first.camera_center_world);
    const auto second_view = Normalize(world_point - second.camera_center_world);
    const auto cosine = std::clamp(Dot(first_view, second_view), -1.0, 1.0);
    const auto angular_separation_degrees = std::acos(cosine) * (180.0 / std::numbers::pi);
    return angular_separation_degrees >= options.min_angular_separation_degrees;
}

bool HasTrackParallax(
    const Vector3d& world_point,
    const std::vector<PointObservationCandidate>& observations,
    const LidarRefinementOptions& options)
{
    for (std::size_t first_index = 0; first_index < observations.size(); ++first_index)
    {
        for (std::size_t second_index = first_index + 1; second_index < observations.size(); ++second_index)
        {
            if (HasGeometricSpread(world_point, observations[first_index], observations[second_index], options))
            {
                return true;
            }
        }
    }

    return false;
}

std::vector<LidarObservationCandidate> BuildImageCandidates(
    const LoadedCapture& capture,
    const std::vector<Vector3d>& world_points,
    const LidarRefinementOptions& options)
{
    if (!capture.depth_info.has_value())
    {
        return {};
    }

    const auto depth_values = LoadDepthMap(*capture.depth_info);
    const auto bucket_width =
        (capture.depth_info->width + options.depth_bucket_size - 1) / options.depth_bucket_size;
    std::unordered_map<int, LidarObservationCandidate> best_by_cell;
    best_by_cell.reserve(
        static_cast<std::size_t>(capture.depth_info->width) *
        static_cast<std::size_t>(capture.depth_info->height) /
        16U);

    for (std::size_t point_index = 0; point_index < world_points.size(); ++point_index)
    {
        const auto projection = ProjectWorldPoint(capture, world_points[point_index]);
        if (!projection.has_value())
        {
            continue;
        }

        if (projection->pixel_x < 0.0 ||
            projection->pixel_x >= static_cast<double>(capture.width) ||
            projection->pixel_y < 0.0 ||
            projection->pixel_y >= static_cast<double>(capture.height))
        {
            continue;
        }

        const auto depth_pixel = MapImagePixelToDepthPixel(capture, projection->pixel_x, projection->pixel_y);
        if (!depth_pixel.has_value())
        {
            continue;
        }

        const auto depth_index = static_cast<std::size_t>(depth_pixel->y) *
                                     static_cast<std::size_t>(capture.depth_info->width) +
                                 static_cast<std::size_t>(depth_pixel->x);
        const auto measured_depth = depth_values[depth_index];
        if (!(measured_depth > 0.0F) || !std::isfinite(measured_depth))
        {
            continue;
        }

        const auto abs_error = std::abs(static_cast<double>(measured_depth) - projection->depth_meters);
        const auto allowed_error = std::max(
            options.base_depth_tolerance_meters,
            static_cast<double>(measured_depth) * options.relative_depth_tolerance);
        if (abs_error > allowed_error)
        {
            continue;
        }

        const auto bucket_x = depth_pixel->x / options.depth_bucket_size;
        const auto bucket_y = depth_pixel->y / options.depth_bucket_size;
        const auto bucket_index = (bucket_y * bucket_width) + bucket_x;
        const auto candidate = LidarObservationCandidate{
            static_cast<int>(point_index),
            projection->pixel_x,
            projection->pixel_y,
            projection->depth_meters,
            abs_error,
            bucket_index};

        const auto it = best_by_cell.find(bucket_index);
        if (it == best_by_cell.end() ||
            candidate.depth_abs_error < it->second.depth_abs_error ||
            (std::abs(candidate.depth_abs_error - it->second.depth_abs_error) < 1e-6 &&
             candidate.depth_meters < it->second.depth_meters))
        {
            best_by_cell[bucket_index] = candidate;
        }
    }

    std::vector<LidarObservationCandidate> candidates;
    candidates.reserve(best_by_cell.size());
    for (const auto& [_, candidate] : best_by_cell)
    {
        candidates.push_back(candidate);
    }

    return candidates;
}

std::vector<PointObservationCandidate> RefineTrack(
    const Vector3d& world_point,
    std::vector<PointObservationCandidate> observations,
    const LidarRefinementOptions& options)
{
    std::sort(
        observations.begin(),
        observations.end(),
        [](const PointObservationCandidate& lhs, const PointObservationCandidate& rhs) {
            if (lhs.depth_abs_error != rhs.depth_abs_error)
            {
                return lhs.depth_abs_error < rhs.depth_abs_error;
            }

            return lhs.capture_index < rhs.capture_index;
        });

    std::vector<PointObservationCandidate> selected;
    selected.reserve(std::min<std::size_t>(observations.size(), static_cast<std::size_t>(options.max_observations_per_track)));
    for (const auto& observation : observations)
    {
        if (selected.empty())
        {
            selected.push_back(observation);
            continue;
        }

        if (selected.size() >= static_cast<std::size_t>(options.max_observations_per_track))
        {
            break;
        }

        const auto has_capture_spread = std::all_of(
            selected.begin(),
            selected.end(),
            [&](const PointObservationCandidate& existing) {
                return std::abs(existing.capture_index - observation.capture_index) >= options.min_capture_gap;
            });
        const auto has_geometric_spread = std::all_of(
            selected.begin(),
            selected.end(),
            [&](const PointObservationCandidate& existing) {
                return HasGeometricSpread(world_point, existing, observation, options);
            });
        if (has_capture_spread || has_geometric_spread)
        {
            selected.push_back(observation);
        }
    }

    if (selected.size() < 2U)
    {
        return {};
    }

    if (!HasTrackParallax(world_point, selected, options))
    {
        return {};
    }

    return selected;
}

ReconstructionBuildResult BuildSparseReconstruction(const LoadedDataset& dataset, const LidarRefinementOptions& options)
{
    std::vector<std::vector<LidarObservationCandidate>> image_candidates;
    image_candidates.reserve(dataset.captures.size());
    std::vector<std::vector<PointObservationCandidate>> point_observations(dataset.world_points.size());

    for (const auto& capture : dataset.captures)
    {
        auto candidates = BuildImageCandidates(capture, dataset.world_points, options);
        for (const auto& candidate : candidates)
        {
            point_observations[static_cast<std::size_t>(candidate.point_index)].push_back(
                PointObservationCandidate{
                    capture.index,
                    capture.index,
                    candidate.pixel_x,
                    candidate.pixel_y,
                    candidate.depth_abs_error,
                    capture.camera_center_world});
        }

        image_candidates.push_back(std::move(candidates));
    }

    std::vector<std::vector<PointObservationCandidate>> refined_observations_by_point(dataset.world_points.size());
    std::vector<std::int64_t> point_id_by_index(dataset.world_points.size(), 0);
    std::int64_t next_point_id = 1;
    for (std::size_t point_index = 0; point_index < point_observations.size(); ++point_index)
    {
        auto& observations = point_observations[point_index];
        if (observations.size() < 2U)
        {
            continue;
        }

        auto refined_track = RefineTrack(dataset.world_points[point_index], observations, options);
        if (refined_track.size() < 2U)
        {
            continue;
        }

        refined_observations_by_point[point_index] = std::move(refined_track);
        point_id_by_index[point_index] = next_point_id++;
    }

    std::vector<std::vector<PointTrackElement>> track_builders(dataset.world_points.size());
    std::vector<double> point_error_sums(dataset.world_points.size(), 0.0);
    std::vector<int> point_error_counts(dataset.world_points.size(), 0);
    std::vector<SparseImagePoints> image_points;
    image_points.reserve(dataset.captures.size());

    for (std::size_t image_index = 0; image_index < dataset.captures.size(); ++image_index)
    {
        const auto& capture = dataset.captures[image_index];
        auto candidates = image_candidates[image_index];
        std::sort(
            candidates.begin(),
            candidates.end(),
            [](const LidarObservationCandidate& lhs, const LidarObservationCandidate& rhs) {
                return lhs.point_index < rhs.point_index;
            });

        std::vector<ImageObservation> observations;
        observations.reserve(candidates.size());
        for (const auto& candidate : candidates)
        {
            const auto point_id = point_id_by_index[static_cast<std::size_t>(candidate.point_index)];
            if (point_id == 0)
            {
                continue;
            }

            const auto& refined_track = refined_observations_by_point[static_cast<std::size_t>(candidate.point_index)];
            const auto keep_for_image = std::any_of(
                refined_track.begin(),
                refined_track.end(),
                [&](const PointObservationCandidate& item) { return item.image_id == capture.index; });
            if (!keep_for_image)
            {
                continue;
            }

            const auto point2d_index = static_cast<int>(observations.size());
            observations.push_back({candidate.pixel_x, candidate.pixel_y, point_id});
            track_builders[static_cast<std::size_t>(candidate.point_index)].push_back({capture.index, point2d_index});
            point_error_sums[static_cast<std::size_t>(candidate.point_index)] += candidate.depth_abs_error;
            point_error_counts[static_cast<std::size_t>(candidate.point_index)] += 1;
        }

        image_points.push_back({capture.index, std::move(observations)});
    }

    std::vector<SparsePoint3D> points3d;
    for (std::size_t point_index = 0; point_index < point_id_by_index.size(); ++point_index)
    {
        const auto point_id = point_id_by_index[point_index];
        if (point_id == 0)
        {
            continue;
        }

        const auto& track = track_builders[point_index];
        if (track.size() < 2U)
        {
            continue;
        }

        const auto error = point_error_counts[point_index] == 0
            ? 0.0
            : point_error_sums[point_index] / static_cast<double>(point_error_counts[point_index]);
        points3d.push_back({
            point_id,
            dataset.world_points[point_index],
            255,
            255,
            255,
            error,
            track});
    }

    return {
        {std::move(image_points), std::move(points3d)},
        dataset.world_points};
}

void WritePointCloudPly(const std::vector<Vector3d>& points, const fs::path& destination_path)
{
    std::ofstream stream(destination_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + destination_path.string() + "' for writing.");
    }

    stream << "ply\n";
    stream << "format ascii 1.0\n";
    stream << "element vertex " << points.size() << "\n";
    stream << "property float x\n";
    stream << "property float y\n";
    stream << "property float z\n";
    stream << "end_header\n";
    for (const auto& point : points)
    {
        stream << ToRoundTrip(point.x) << ' ' << ToRoundTrip(point.y) << ' ' << ToRoundTrip(point.z) << '\n';
    }
}

void WriteCamerasText(const std::vector<LoadedCapture>& captures, const fs::path& destination_path)
{
    std::ofstream stream(destination_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + destination_path.string() + "' for writing.");
    }

    stream << "# Camera list with one line of data per camera:\n";
    stream << "#   CAMERA_ID, MODEL, WIDTH, HEIGHT, PARAMS[]\n";
    stream << "# Number of cameras: " << captures.size() << "\n";
    for (const auto& capture : captures)
    {
        stream << capture.index << " PINHOLE "
               << capture.width << ' '
               << capture.height << ' '
               << ToRoundTrip(capture.fx) << ' '
               << ToRoundTrip(capture.fy) << ' '
               << ToRoundTrip(capture.cx) << ' '
               << ToRoundTrip(capture.cy) << '\n';
    }
}

void WriteImagesText(
    const std::vector<LoadedCapture>& captures,
    const SparseReconstruction& sparse_reconstruction,
    const fs::path& destination_path)
{
    std::map<int, const SparseImagePoints*> observations_by_image;
    for (const auto& image : sparse_reconstruction.image_points)
    {
        observations_by_image[image.image_id] = &image;
    }

    std::ofstream stream(destination_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + destination_path.string() + "' for writing.");
    }

    stream << "# Image list with two lines of data per image:\n";
    stream << "#   IMAGE_ID, QW, QX, QY, QZ, TX, TY, TZ, CAMERA_ID, IMAGE_NAME\n";
    stream << "#   POINTS2D[] as (X, Y, POINT3D_ID)\n";
    stream << "# Number of images: " << captures.size() << "\n";

    for (const auto& capture : captures)
    {
        const auto& q = capture.colmap_pose.rotation;
        const auto& t = capture.colmap_pose.translation;
        stream << capture.index << ' '
               << ToRoundTrip(q.w) << ' '
               << ToRoundTrip(q.x) << ' '
               << ToRoundTrip(q.y) << ' '
               << ToRoundTrip(q.z) << ' '
               << ToRoundTrip(t.x) << ' '
               << ToRoundTrip(t.y) << ' '
               << ToRoundTrip(t.z) << ' '
               << capture.index << ' '
               << capture.image_name << '\n';

        const auto it = observations_by_image.find(capture.index);
        if (it == observations_by_image.end() || it->second->observations.empty())
        {
            stream << '\n';
            continue;
        }

        for (std::size_t observation_index = 0; observation_index < it->second->observations.size(); ++observation_index)
        {
            const auto& observation = it->second->observations[observation_index];
            if (observation_index > 0U)
            {
                stream << ' ';
            }

            stream << ToRoundTrip(observation.x) << ' '
                   << ToRoundTrip(observation.y) << ' '
                   << observation.point3d_id;
        }
        stream << '\n';
    }
}

void WritePoints3DText(const SparseReconstruction& sparse_reconstruction, const fs::path& destination_path)
{
    std::ofstream stream(destination_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + destination_path.string() + "' for writing.");
    }

    stream << "# 3D point list with one line of data per point:\n";
    stream << "#   POINT3D_ID, X, Y, Z, R, G, B, ERROR, TRACK[] as (IMAGE_ID, POINT2D_IDX)\n";
    stream << "# Number of points: " << sparse_reconstruction.points3d.size() << "\n";
    for (const auto& point : sparse_reconstruction.points3d)
    {
        stream << point.point3d_id << ' '
               << ToRoundTrip(point.position.x) << ' '
               << ToRoundTrip(point.position.y) << ' '
               << ToRoundTrip(point.position.z) << ' '
               << static_cast<int>(point.r) << ' '
               << static_cast<int>(point.g) << ' '
               << static_cast<int>(point.b) << ' '
               << ToRoundTrip(point.error);
        for (const auto& track_element : point.track)
        {
            stream << ' ' << track_element.image_id << ' ' << track_element.point2d_index;
        }
        stream << '\n';
    }
}

void WriteCamerasBinary(const std::vector<LoadedCapture>& captures, const fs::path& destination_path)
{
    std::ofstream stream(destination_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + destination_path.string() + "' for writing.");
    }

    WriteLittleEndian<std::uint64_t>(stream, static_cast<std::uint64_t>(captures.size()));
    for (const auto& capture : captures)
    {
        WriteLittleEndian<std::int32_t>(stream, capture.index);
        WriteLittleEndian<std::int32_t>(stream, 1);
        WriteLittleEndian<std::uint64_t>(stream, static_cast<std::uint64_t>(capture.width));
        WriteLittleEndian<std::uint64_t>(stream, static_cast<std::uint64_t>(capture.height));
        WriteLittleEndian<double>(stream, capture.fx);
        WriteLittleEndian<double>(stream, capture.fy);
        WriteLittleEndian<double>(stream, capture.cx);
        WriteLittleEndian<double>(stream, capture.cy);
    }
}

void WriteImagesBinary(
    const std::vector<LoadedCapture>& captures,
    const SparseReconstruction& sparse_reconstruction,
    const fs::path& destination_path)
{
    std::map<int, const SparseImagePoints*> observations_by_image;
    for (const auto& image : sparse_reconstruction.image_points)
    {
        observations_by_image[image.image_id] = &image;
    }

    std::ofstream stream(destination_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + destination_path.string() + "' for writing.");
    }

    WriteLittleEndian<std::uint64_t>(stream, static_cast<std::uint64_t>(captures.size()));
    for (const auto& capture : captures)
    {
        const auto& q = capture.colmap_pose.rotation;
        const auto& t = capture.colmap_pose.translation;
        WriteLittleEndian<std::int32_t>(stream, capture.index);
        WriteLittleEndian<double>(stream, q.w);
        WriteLittleEndian<double>(stream, q.x);
        WriteLittleEndian<double>(stream, q.y);
        WriteLittleEndian<double>(stream, q.z);
        WriteLittleEndian<double>(stream, t.x);
        WriteLittleEndian<double>(stream, t.y);
        WriteLittleEndian<double>(stream, t.z);
        WriteLittleEndian<std::int32_t>(stream, capture.index);
        stream.write(capture.image_name.data(), static_cast<std::streamsize>(capture.image_name.size()));
        stream.put('\0');

        const auto it = observations_by_image.find(capture.index);
        if (it == observations_by_image.end())
        {
            WriteLittleEndian<std::uint64_t>(stream, 0);
            continue;
        }

        WriteLittleEndian<std::uint64_t>(stream, static_cast<std::uint64_t>(it->second->observations.size()));
        for (const auto& observation : it->second->observations)
        {
            WriteLittleEndian<double>(stream, observation.x);
            WriteLittleEndian<double>(stream, observation.y);
            WriteLittleEndian<std::int64_t>(stream, observation.point3d_id);
        }
    }
}

void WritePoints3DBinary(const SparseReconstruction& sparse_reconstruction, const fs::path& destination_path)
{
    std::ofstream stream(destination_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + destination_path.string() + "' for writing.");
    }

    WriteLittleEndian<std::uint64_t>(stream, static_cast<std::uint64_t>(sparse_reconstruction.points3d.size()));
    for (const auto& point : sparse_reconstruction.points3d)
    {
        WriteLittleEndian<std::uint64_t>(stream, static_cast<std::uint64_t>(point.point3d_id));
        WriteLittleEndian<double>(stream, point.position.x);
        WriteLittleEndian<double>(stream, point.position.y);
        WriteLittleEndian<double>(stream, point.position.z);
        WriteLittleEndian<std::uint8_t>(stream, point.r);
        WriteLittleEndian<std::uint8_t>(stream, point.g);
        WriteLittleEndian<std::uint8_t>(stream, point.b);
        WriteLittleEndian<double>(stream, point.error);
        WriteLittleEndian<std::uint64_t>(stream, static_cast<std::uint64_t>(point.track.size()));
        for (const auto& track_element : point.track)
        {
            WriteLittleEndian<std::int32_t>(stream, track_element.image_id);
            WriteLittleEndian<std::int32_t>(stream, track_element.point2d_index);
        }
    }
}

void CopyImages(const std::vector<LoadedCapture>& captures, const fs::path& images_directory)
{
    for (const auto& capture : captures)
    {
        const auto destination_path = images_directory / capture.image_name;
        fs::copy_file(capture.image_path, destination_path, fs::copy_options::overwrite_existing);
    }
}

void WriteColmapOutput(
    const LoadedDataset& dataset,
    const ReconstructionBuildResult& reconstruction,
    const fs::path& output_directory,
    const bool copy_images)
{
    fs::create_directories(output_directory);
    const auto sparse_directory = output_directory / "sparse" / "0";
    const auto images_directory = output_directory / "images";
    fs::create_directories(sparse_directory);
    fs::create_directories(images_directory);

    const auto& sparse_reconstruction = reconstruction.sparse_reconstruction;
    WriteCamerasText(dataset.captures, sparse_directory / "cameras.txt");
    WriteImagesText(dataset.captures, sparse_reconstruction, sparse_directory / "images.txt");
    WritePoints3DText(sparse_reconstruction, sparse_directory / "points3D.txt");
    WriteCamerasBinary(dataset.captures, sparse_directory / "cameras.bin");
    WriteImagesBinary(dataset.captures, sparse_reconstruction, sparse_directory / "images.bin");
    WritePoints3DBinary(sparse_reconstruction, sparse_directory / "points3D.bin");
    WritePointCloudPly(dataset.world_points, output_directory / "lidar_dense_world_points.ply");

    std::vector<Vector3d> sparse_points;
    sparse_points.reserve(sparse_reconstruction.points3d.size());
    for (const auto& point : sparse_reconstruction.points3d)
    {
        sparse_points.push_back(point.position);
    }
    WritePointCloudPly(sparse_points, output_directory / "lidar_sparse_world_points.ply");

    if (copy_images)
    {
        CopyImages(dataset.captures, images_directory);
    }
}

std::uint64_t ReadTextCameraCount(const fs::path& path)
{
    std::ifstream stream(path);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + path.string() + "'.");
    }

    std::uint64_t count = 0;
    std::string line;
    while (std::getline(stream, line))
    {
        if (line.empty() || (!line.empty() && line[0] == '#'))
        {
            continue;
        }

        ++count;
    }

    return count;
}

std::pair<std::uint64_t, std::uint64_t> ReadTextImageCounts(const fs::path& path)
{
    std::ifstream stream(path);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + path.string() + "'.");
    }

    std::uint64_t image_count = 0;
    std::uint64_t observation_count = 0;
    auto expecting_header = true;
    std::string line;
    while (std::getline(stream, line))
    {
        if (line.empty() || (!line.empty() && line[0] == '#'))
        {
            continue;
        }

        const auto parts = SplitWhitespace(line);
        if (expecting_header)
        {
            ++image_count;
        }
        else
        {
            observation_count += static_cast<std::uint64_t>(parts.size() / 3U);
        }

        expecting_header = !expecting_header;
    }

    return {image_count, observation_count};
}

std::pair<std::uint64_t, std::uint64_t> ReadTextPointCounts(const fs::path& path)
{
    std::ifstream stream(path);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + path.string() + "'.");
    }

    std::uint64_t point_count = 0;
    std::uint64_t track_element_count = 0;
    std::string line;
    while (std::getline(stream, line))
    {
        if (line.empty() || (!line.empty() && line[0] == '#'))
        {
            continue;
        }

        ++point_count;
        const auto parts = SplitWhitespace(line);
        if (parts.size() >= 8U)
        {
            track_element_count += static_cast<std::uint64_t>((parts.size() - 8U) / 2U);
        }
    }

    return {point_count, track_element_count};
}

std::string ReadNullTerminatedString(std::istream& stream)
{
    std::string value;
    for (;;)
    {
        char next = 0;
        stream.read(&next, 1);
        if (!stream)
        {
            throw std::runtime_error("Unexpected end of stream while reading a string.");
        }

        if (next == '\0')
        {
            break;
        }

        value.push_back(next);
    }

    return value;
}

std::pair<std::uint64_t, std::uint64_t> ReadBinaryImageCounts(const fs::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + path.string() + "'.");
    }

    const auto image_count = ReadLittleEndian<std::uint64_t>(stream);
    std::uint64_t observation_count = 0;
    for (std::uint64_t image_index = 0; image_index < image_count; ++image_index)
    {
        static_cast<void>(ReadLittleEndian<std::int32_t>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<std::int32_t>(stream));
        static_cast<void>(ReadNullTerminatedString(stream));

        const auto points2d = ReadLittleEndian<std::uint64_t>(stream);
        observation_count += points2d;
        for (std::uint64_t point_index = 0; point_index < points2d; ++point_index)
        {
            static_cast<void>(ReadLittleEndian<double>(stream));
            static_cast<void>(ReadLittleEndian<double>(stream));
            static_cast<void>(ReadLittleEndian<std::int64_t>(stream));
        }
    }

    return {image_count, observation_count};
}

std::pair<std::uint64_t, std::uint64_t> ReadBinaryPointCounts(const fs::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + path.string() + "'.");
    }

    const auto point_count = ReadLittleEndian<std::uint64_t>(stream);
    std::uint64_t track_element_count = 0;
    for (std::uint64_t point_index = 0; point_index < point_count; ++point_index)
    {
        static_cast<void>(ReadLittleEndian<std::uint64_t>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        static_cast<void>(ReadLittleEndian<std::uint8_t>(stream));
        static_cast<void>(ReadLittleEndian<std::uint8_t>(stream));
        static_cast<void>(ReadLittleEndian<std::uint8_t>(stream));
        static_cast<void>(ReadLittleEndian<double>(stream));
        const auto track_length = ReadLittleEndian<std::uint64_t>(stream);
        track_element_count += track_length;
        for (std::uint64_t track_index = 0; track_index < track_length; ++track_index)
        {
            static_cast<void>(ReadLittleEndian<std::int32_t>(stream));
            static_cast<void>(ReadLittleEndian<std::int32_t>(stream));
        }
    }

    return {point_count, track_element_count};
}

void WriteModelValidationSummary(const fs::path& output_directory, const ColmapModelValidationSummary& summary)
{
    std::ofstream stream(output_directory / "model_validation.txt", std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to write model_validation.txt.");
    }

    stream << "COLMAP export validation\n";
    stream << "match=" << BoolText(summary.IsMatch()) << "\n";
    stream << "text_cameras=" << summary.text_counts.cameras << "\n";
    stream << "bin_cameras=" << summary.binary_counts.cameras << "\n";
    stream << "text_images=" << summary.text_counts.images << "\n";
    stream << "bin_images=" << summary.binary_counts.images << "\n";
    stream << "text_points3D=" << summary.text_counts.points3d << "\n";
    stream << "bin_points3D=" << summary.binary_counts.points3d << "\n";
    stream << "text_image_observations=" << summary.text_counts.image_observations << "\n";
    stream << "bin_image_observations=" << summary.binary_counts.image_observations << "\n";
    stream << "text_point_track_elements=" << summary.text_counts.point_track_elements << "\n";
    stream << "bin_point_track_elements=" << summary.binary_counts.point_track_elements << "\n";
}

ColmapModelValidationSummary ValidateExportedModel(const fs::path& output_directory)
{
    const auto sparse_directory = output_directory / "sparse" / "0";
    const auto text_counts = ColmapModelCounts{
        ReadTextCameraCount(sparse_directory / "cameras.txt"),
        ReadTextImageCounts(sparse_directory / "images.txt").first,
        ReadTextPointCounts(sparse_directory / "points3D.txt").first,
        ReadTextImageCounts(sparse_directory / "images.txt").second,
        ReadTextPointCounts(sparse_directory / "points3D.txt").second};
    const auto binary_counts = ColmapModelCounts{
        [&]() {
            std::ifstream stream(sparse_directory / "cameras.bin", std::ios::binary);
            if (!stream)
            {
                throw std::runtime_error("Failed to open cameras.bin.");
            }
            return ReadLittleEndian<std::uint64_t>(stream);
        }(),
        ReadBinaryImageCounts(sparse_directory / "images.bin").first,
        ReadBinaryPointCounts(sparse_directory / "points3D.bin").first,
        ReadBinaryImageCounts(sparse_directory / "images.bin").second,
        ReadBinaryPointCounts(sparse_directory / "points3D.bin").second};

    const auto summary = ColmapModelValidationSummary{text_counts, binary_counts};
    WriteModelValidationSummary(output_directory, summary);
    return summary;
}

Vector3d Average(const std::vector<Vector3d>& points)
{
    Vector3d sum{};
    for (const auto& point : points)
    {
        sum += point;
    }

    return points.empty() ? Vector3d{} : (sum / static_cast<double>(points.size()));
}

SimilarityTransform EstimateSimilarity(const std::vector<std::pair<Vector3d, Vector3d>>& matches)
{
    std::vector<Vector3d> source_points;
    std::vector<Vector3d> target_points;
    source_points.reserve(matches.size());
    target_points.reserve(matches.size());
    for (const auto& match : matches)
    {
        source_points.push_back(match.first);
        target_points.push_back(match.second);
    }

    const auto source_centroid = Average(source_points);
    const auto target_centroid = Average(target_points);

    std::array<std::array<double, 3>, 3> covariance{};
    for (const auto& match : matches)
    {
        const auto source = match.first - source_centroid;
        const auto target = match.second - target_centroid;
        covariance[0][0] += source.x * target.x;
        covariance[0][1] += source.x * target.y;
        covariance[0][2] += source.x * target.z;
        covariance[1][0] += source.y * target.x;
        covariance[1][1] += source.y * target.y;
        covariance[1][2] += source.y * target.z;
        covariance[2][0] += source.z * target.x;
        covariance[2][1] += source.z * target.y;
        covariance[2][2] += source.z * target.z;
    }

    const auto sxx = covariance[0][0];
    const auto sxy = covariance[0][1];
    const auto sxz = covariance[0][2];
    const auto syx = covariance[1][0];
    const auto syy = covariance[1][1];
    const auto syz = covariance[1][2];
    const auto szx = covariance[2][0];
    const auto szy = covariance[2][1];
    const auto szz = covariance[2][2];

    std::array<std::array<double, 4>, 4> horn{};
    horn[0][0] = sxx + syy + szz;
    horn[0][1] = syz - szy;
    horn[0][2] = szx - sxz;
    horn[0][3] = sxy - syx;
    horn[1][0] = syz - szy;
    horn[1][1] = sxx - syy - szz;
    horn[1][2] = sxy + syx;
    horn[1][3] = szx + sxz;
    horn[2][0] = szx - sxz;
    horn[2][1] = sxy + syx;
    horn[2][2] = -sxx + syy - szz;
    horn[2][3] = syz + szy;
    horn[3][0] = sxy - syx;
    horn[3][1] = szx + sxz;
    horn[3][2] = syz + szy;
    horn[3][3] = -sxx - syy + szz;

    std::array<double, 4> eigenvector{1.0, 0.0, 0.0, 0.0};
    for (int iteration = 0; iteration < 64; ++iteration)
    {
        std::array<double, 4> next{};
        for (int row = 0; row < 4; ++row)
        {
            for (int col = 0; col < 4; ++col)
            {
                next[row] += horn[row][col] * eigenvector[col];
            }
        }

        const auto norm =
            std::sqrt((next[0] * next[0]) + (next[1] * next[1]) + (next[2] * next[2]) + (next[3] * next[3]));
        if (norm <= 1e-12)
        {
            break;
        }

        for (int index = 0; index < 4; ++index)
        {
            eigenvector[index] = next[index] / norm;
        }
    }

    auto quaternion = Normalize(Quaterniond{eigenvector[0], eigenvector[1], eigenvector[2], eigenvector[3]});
    if (quaternion.w < 0.0)
    {
        quaternion = Negate(quaternion);
    }

    const auto rotation = QuaternionToRotationMatrix(quaternion);
    auto numerator = 0.0;
    auto denominator = 0.0;
    for (const auto& match : matches)
    {
        const auto source = match.first - source_centroid;
        const auto target = match.second - target_centroid;
        const auto rotated = Multiply(rotation, source);
        numerator += Dot(target, rotated);
        denominator += Dot(source, source);
    }

    const auto scale = denominator <= 1e-12 ? 1.0 : (numerator / denominator);
    const auto translation = target_centroid - (scale * Multiply(rotation, source_centroid));
    return {rotation, translation, scale};
}

Vector3d TransformPoint(const SimilarityTransform& transform, const Vector3d& point)
{
    return (transform.scale * Multiply(transform.rotation, point)) + transform.translation;
}

ColmapPose TransformPose(const ColmapPose& pose, const Vector3d& original_center, const SimilarityTransform& transform)
{
    const auto transformed_center = TransformPoint(transform, original_center);
    const auto transformed_rotation = Multiply(pose.world_to_camera, Transpose(transform.rotation));
    const auto transformed_translation = Vector3d{
        -Dot({transformed_rotation.m11, transformed_rotation.m12, transformed_rotation.m13}, transformed_center),
        -Dot({transformed_rotation.m21, transformed_rotation.m22, transformed_rotation.m23}, transformed_center),
        -Dot({transformed_rotation.m31, transformed_rotation.m32, transformed_rotation.m33}, transformed_center)};

    auto quaternion = RotationMatrixToHamiltonQuaternion(transformed_rotation);
    if (quaternion.w < 0.0)
    {
        quaternion = Negate(quaternion);
    }

    return {transformed_rotation, quaternion, transformed_translation};
}

ColmapPose DeriveMetadataPose(const ColmapPose& display_pose, const ImageRotation image_rotation)
{
    switch (image_rotation)
    {
        case ImageRotation::Rotate90Clockwise:
            return RotateImage90CounterClockwise(display_pose);
        case ImageRotation::Rotate90CounterClockwise:
            return RotateImage90Clockwise(display_pose);
        case ImageRotation::None:
        default:
            return display_pose;
    }
}

std::string ResolveModelDirectory(const fs::path& reference_model_path)
{
    const auto full_path = fs::absolute(reference_model_path);
    if (fs::is_directory(full_path) && fs::exists(full_path / "images.txt"))
    {
        return full_path.string();
    }

    const auto sparse_candidate = full_path / "sparse" / "0";
    if (fs::is_directory(sparse_candidate) && fs::exists(sparse_candidate / "images.txt"))
    {
        return sparse_candidate.string();
    }

    throw std::runtime_error("Could not locate a COLMAP text model under '" + reference_model_path.string() + "'.");
}

std::unordered_map<std::string, ReferenceImagePose> LoadReferencePoses(const fs::path& images_path)
{
    std::ifstream stream(images_path);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + images_path.string() + "'.");
    }

    std::unordered_map<std::string, ReferenceImagePose> poses;
    auto expecting_image_header = true;
    std::string line;
    while (std::getline(stream, line))
    {
        if (!line.empty() && line[0] == '#')
        {
            continue;
        }

        if (!expecting_image_header)
        {
            expecting_image_header = true;
            continue;
        }

        if (line.empty())
        {
            continue;
        }

        const auto parts = SplitWhitespace(line);
        if (parts.size() >= 10U)
        {
            const auto qw = std::stod(parts[1]);
            const auto qx = std::stod(parts[2]);
            const auto qy = std::stod(parts[3]);
            const auto qz = std::stod(parts[4]);
            const auto tx = std::stod(parts[5]);
            const auto ty = std::stod(parts[6]);
            const auto tz = std::stod(parts[7]);
            const auto image_name = parts[9];

            auto rotation_quaternion = Normalize(Quaterniond{qw, qx, qy, qz});
            auto rotation = QuaternionToRotationMatrix(rotation_quaternion);
            const auto translation = Vector3d{tx, ty, tz};
            const auto center = Vector3d{
                -((rotation.m11 * translation.x) + (rotation.m21 * translation.y) + (rotation.m31 * translation.z)),
                -((rotation.m12 * translation.x) + (rotation.m22 * translation.y) + (rotation.m32 * translation.z)),
                -((rotation.m13 * translation.x) + (rotation.m23 * translation.y) + (rotation.m33 * translation.z))};

            poses[image_name] = {image_name, center, {rotation, rotation_quaternion, translation}};
        }

        expecting_image_header = false;
    }

    return poses;
}

std::tuple<std::int64_t, int, std::streamoff> ReadBinaryPlyHeader(const fs::path& ply_path)
{
    std::ifstream stream(ply_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + ply_path.string() + "'.");
    }

    std::string current_line;
    std::int64_t vertex_count = 0;
    auto property_count = 0;
    auto in_vertex_element = false;
    while (true)
    {
        char next_byte = 0;
        stream.read(&next_byte, 1);
        if (!stream)
        {
            throw std::runtime_error("Unexpected end of file while reading PLY header '" + ply_path.string() + "'.");
        }

        if (next_byte == '\n')
        {
            if (!current_line.empty() && current_line.back() == '\r')
            {
                current_line.pop_back();
            }

            if (current_line.rfind("format ", 0) == 0)
            {
                if (current_line != "format binary_little_endian 1.0")
                {
                    throw std::runtime_error("Unsupported PLY format in '" + ply_path.string() + "': " + current_line);
                }
            }
            else if (current_line.rfind("element ", 0) == 0)
            {
                const auto parts = SplitWhitespace(current_line);
                in_vertex_element = parts.size() >= 3U && parts[1] == "vertex";
                if (in_vertex_element)
                {
                    vertex_count = std::stoll(parts[2]);
                }
            }
            else if (current_line.rfind("property ", 0) == 0 && in_vertex_element)
            {
                ++property_count;
            }
            else if (current_line == "end_header")
            {
                return {vertex_count, property_count, stream.tellg()};
            }

            current_line.clear();
            continue;
        }

        current_line.push_back(next_byte);
    }
}

std::vector<Vector3d> LoadBinaryPlySample(const fs::path& ply_path, const int max_sample_count)
{
    auto [vertex_count, property_count, data_start_offset] = ReadBinaryPlyHeader(ply_path);
    if (vertex_count <= 0 || property_count < 3)
    {
        return {};
    }

    std::ifstream stream(ply_path, std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to open '" + ply_path.string() + "'.");
    }

    stream.seekg(data_start_offset);
    const auto stride = std::max(1, static_cast<int>(std::ceil(static_cast<double>(vertex_count) / static_cast<double>(max_sample_count))));
    const auto record_size_bytes = property_count * static_cast<int>(sizeof(float));
    const auto trailing_bytes = record_size_bytes - static_cast<int>(3 * sizeof(float));
    std::vector<Vector3d> points;
    points.reserve(static_cast<std::size_t>(std::min<std::int64_t>(vertex_count / stride, max_sample_count)));

    for (std::int64_t vertex_index = 0; vertex_index < vertex_count; ++vertex_index)
    {
        const auto x = ReadLittleEndian<float>(stream);
        const auto y = ReadLittleEndian<float>(stream);
        const auto z = ReadLittleEndian<float>(stream);
        if (trailing_bytes > 0)
        {
            stream.seekg(trailing_bytes, std::ios::cur);
        }

        if ((vertex_index % stride) != 0)
        {
            continue;
        }

        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        {
            continue;
        }

        points.push_back({static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)});
    }

    return points;
}

std::vector<Vector3d> LoadReferencePoints(const fs::path& resolved_model_directory)
{
    const auto root_directory = fs::absolute(resolved_model_directory / ".." / "..");
    const auto splat_ply_path = root_directory / "output" / "splat_30000.ply";
    if (fs::exists(splat_ply_path))
    {
        auto dense_points = LoadBinaryPlySample(splat_ply_path, 200000);
        if (!dense_points.empty())
        {
            return dense_points;
        }
    }

    const auto sparse_points_path = resolved_model_directory / "points3D.txt";
    std::ifstream stream(sparse_points_path);
    if (!stream)
    {
        throw std::runtime_error("Could not load reference points from '" + sparse_points_path.string() + "'.");
    }

    std::vector<Vector3d> sparse_points;
    std::string line;
    while (std::getline(stream, line))
    {
        if (line.empty() || (!line.empty() && line[0] == '#'))
        {
            continue;
        }

        const auto parts = SplitWhitespace(line);
        if (parts.size() < 4U)
        {
            continue;
        }

        sparse_points.push_back({
            std::stod(parts[1]),
            std::stod(parts[2]),
            std::stod(parts[3])});
    }

    if (!sparse_points.empty())
    {
        return sparse_points;
    }

    throw std::runtime_error("Could not load reference points from '" + sparse_points_path.string() + "'.");
}

SimilarityTransform RefineSimilarityWithIcp(
    const std::vector<Vector3d>& source_points,
    const std::vector<Vector3d>& target_points,
    const SimilarityTransform& initial_transform)
{
    if (source_points.empty() || target_points.empty())
    {
        return initial_transform;
    }

    auto current = initial_transform;
    const auto target_index = TargetPointIndex::Build(target_points, 1.0);
    const std::array<double, 5> thresholds{5.0, 3.0, 2.0, 1.0, 0.5};
    constexpr int source_stride = 4;

    for (const auto threshold : thresholds)
    {
        std::vector<std::pair<Vector3d, Vector3d>> correspondences;
        correspondences.reserve(source_points.size() / source_stride);
        for (std::size_t point_index = 0; point_index < source_points.size(); point_index += source_stride)
        {
            const auto transformed = TransformPoint(current, source_points[point_index]);
            Vector3d nearest{};
            if (target_index.TryFindNearest(transformed, threshold, nearest))
            {
                correspondences.emplace_back(source_points[point_index], nearest);
            }
        }

        if (correspondences.size() < 64U)
        {
            continue;
        }

        current = EstimateSimilarity(correspondences);
    }

    return current;
}

AlignmentSummary BuildAlignmentSummary(
    const std::vector<std::tuple<LoadedCapture, Vector3d, ReferenceImagePose>>& matches,
    const std::string& reference_model_path,
    const SimilarityTransform& transform,
    const bool snap_shared_camera_poses_to_reference)
{
    auto sum = 0.0;
    auto max = 0.0;
    for (const auto& match : matches)
    {
        const auto transformed_center = TransformPoint(transform, std::get<1>(match));
        const auto residual = Distance(transformed_center, std::get<2>(match).camera_center_world);
        sum += residual;
        max = std::max(max, residual);
    }

    return {
        static_cast<int>(matches.size()),
        snap_shared_camera_poses_to_reference ? static_cast<int>(matches.size()) : 0,
        snap_shared_camera_poses_to_reference,
        transform.scale,
        matches.empty() ? 0.0 : (sum / static_cast<double>(matches.size())),
        max,
        reference_model_path};
}

bool ShouldUseRefinedTransform(const AlignmentSummary& initial, const AlignmentSummary& refined)
{
    if (!std::isfinite(refined.scale) || refined.scale <= 1e-9)
    {
        return false;
    }

    if (!std::isfinite(refined.mean_center_residual_meters) ||
        !std::isfinite(refined.max_center_residual_meters))
    {
        return false;
    }

    return refined.mean_center_residual_meters < initial.mean_center_residual_meters;
}

LoadedCapture CreateAlignedCapture(
    const LoadedCapture& capture,
    const SimilarityTransform& transform)
{
    auto aligned = capture;
    aligned.camera_center_world = TransformPoint(transform, capture.camera_center_world);
    aligned.metadata_colmap_pose = TransformPose(capture.metadata_colmap_pose, capture.camera_center_world, transform);
    aligned.colmap_pose = TransformPose(capture.colmap_pose, capture.camera_center_world, transform);
    return aligned;
}

LoadedCapture CreateAlignedCapture(
    const LoadedCapture& capture,
    const std::unordered_map<std::string, ReferenceImagePose>& reference_poses,
    const SimilarityTransform& transform,
    const bool snap_shared_camera_poses_to_reference)
{
    if (snap_shared_camera_poses_to_reference)
    {
        const auto it = reference_poses.find(capture.image_name);
        if (it != reference_poses.end())
        {
            auto aligned = capture;
            aligned.camera_center_world = it->second.camera_center_world;
            aligned.metadata_colmap_pose = DeriveMetadataPose(it->second.pose, capture.image_rotation);
            aligned.colmap_pose = it->second.pose;
            return aligned;
        }
    }

    return CreateAlignedCapture(capture, transform);
}

ReconstructionBuildResult TransformReconstruction(
    const ReconstructionBuildResult& reconstruction,
    const SimilarityTransform& transform)
{
    auto transformed = reconstruction;
    for (auto& point : transformed.sparse_reconstruction.points3d)
    {
        point.position = TransformPoint(transform, point.position);
    }

    for (auto& point : transformed.dense_aligned_points)
    {
        point = TransformPoint(transform, point);
    }

    return transformed;
}

void WriteAlignmentSummary(const AlignmentSummary& summary, const fs::path& output_directory)
{
    std::ofstream stream(output_directory / "alignment_summary.txt", std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to write alignment_summary.txt.");
    }

    stream << "Alignment summary\n";
    stream << "reference_model_path=" << summary.reference_model_path << "\n";
    stream << "matched_image_count=" << summary.matched_image_count << "\n";
    stream << "snapped_image_count=" << summary.snapped_image_count << "\n";
    stream << "shared_camera_poses_snapped_to_reference="
           << BoolText(summary.shared_camera_poses_snapped_to_reference) << "\n";
    stream << "scale=" << ToRoundTrip(summary.scale) << "\n";
    stream << "mean_center_residual_m=" << ToRoundTrip(summary.mean_center_residual_meters) << "\n";
    stream << "max_center_residual_m=" << ToRoundTrip(summary.max_center_residual_meters) << "\n";
}

AlignmentResult ApplyReferenceAlignment(
    const LoadedDataset& dataset,
    const ReconstructionBuildResult& reconstruction,
    const fs::path& reference_model_path,
    const bool snap_shared_camera_poses_to_reference)
{
    const auto resolved_model_directory = ResolveModelDirectory(reference_model_path);
    const auto reference_poses = LoadReferencePoses(fs::path(resolved_model_directory) / "images.txt");
    const auto reference_points = LoadReferencePoints(fs::path(resolved_model_directory));

    std::vector<std::tuple<LoadedCapture, Vector3d, ReferenceImagePose>> camera_matches;
    for (const auto& capture : dataset.captures)
    {
        const auto it = reference_poses.find(capture.image_name);
        if (it == reference_poses.end())
        {
            continue;
        }

        camera_matches.emplace_back(capture, capture.camera_center_world, it->second);
    }

    if (camera_matches.size() < 3U)
    {
        throw std::runtime_error(
            "Need at least 3 shared images to align against '" +
            reference_model_path.string() +
            "', found " +
            std::to_string(camera_matches.size()) +
            ".");
    }

    std::vector<std::pair<Vector3d, Vector3d>> center_matches;
    center_matches.reserve(camera_matches.size());
    for (const auto& match : camera_matches)
    {
        center_matches.emplace_back(std::get<1>(match), std::get<2>(match).camera_center_world);
    }

    const auto initial_transform = EstimateSimilarity(center_matches);
    const auto refined_transform = RefineSimilarityWithIcp(dataset.world_points, reference_points, initial_transform);
    const auto initial_summary = BuildAlignmentSummary(
        camera_matches,
        reference_model_path.string(),
        initial_transform,
        false);
    const auto refined_summary = BuildAlignmentSummary(
        camera_matches,
        reference_model_path.string(),
        refined_transform,
        false);
    const auto chosen_transform =
        ShouldUseRefinedTransform(initial_summary, refined_summary) ? refined_transform : initial_transform;

    LoadedDataset aligned_dataset = dataset;
    aligned_dataset.captures.clear();
    aligned_dataset.captures.reserve(dataset.captures.size());
    for (const auto& capture : dataset.captures)
    {
        aligned_dataset.captures.push_back(CreateAlignedCapture(
            capture,
            reference_poses,
            chosen_transform,
            snap_shared_camera_poses_to_reference));
    }

    aligned_dataset.world_points.clear();
    aligned_dataset.world_points.reserve(dataset.world_points.size());
    for (const auto& point : dataset.world_points)
    {
        aligned_dataset.world_points.push_back(TransformPoint(chosen_transform, point));
    }

    const auto aligned_reconstruction = TransformReconstruction(reconstruction, chosen_transform);
    const auto summary = BuildAlignmentSummary(
        camera_matches,
        reference_model_path.string(),
        chosen_transform,
        snap_shared_camera_poses_to_reference);
    return {std::move(aligned_dataset), aligned_reconstruction, summary};
}

void WriteCameraPoseSummary(const CameraPoseSummary& summary, const fs::path& output_directory)
{
    std::ofstream stream(output_directory / "camera_pose_summary.txt", std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to write camera_pose_summary.txt.");
    }

    stream << "Camera pose summary\n";
    stream << "camera_count=" << summary.camera_count << "\n";
    stream << "mean_center_error_m=" << ToRoundTrip(summary.mean_center_error_meters) << "\n";
    stream << "max_center_error_m=" << ToRoundTrip(summary.max_center_error_meters) << "\n";
}

CameraPoseSummary RunCameraPoseValidation(const LoadedDataset& dataset, const fs::path& output_directory)
{
    auto sum = 0.0;
    auto max = 0.0;
    for (const auto& capture : dataset.captures)
    {
        const auto& rotation = capture.colmap_pose.world_to_camera;
        const auto& translation = capture.colmap_pose.translation;
        const auto center = Vector3d{
            -((rotation.m11 * translation.x) + (rotation.m21 * translation.y) + (rotation.m31 * translation.z)),
            -((rotation.m12 * translation.x) + (rotation.m22 * translation.y) + (rotation.m32 * translation.z)),
            -((rotation.m13 * translation.x) + (rotation.m23 * translation.y) + (rotation.m33 * translation.z))};
        const auto error = Distance(center, capture.camera_center_world);
        sum += error;
        max = std::max(max, error);
    }

    const auto summary = CameraPoseSummary{
        static_cast<int>(dataset.captures.size()),
        dataset.captures.empty() ? 0.0 : (sum / static_cast<double>(dataset.captures.size())),
        max};
    WriteCameraPoseSummary(summary, output_directory);
    return summary;
}

std::vector<LoadedCapture> SelectValidationCaptures(const std::vector<LoadedCapture>& captures, const int requested_count)
{
    if (static_cast<int>(captures.size()) <= requested_count)
    {
        return captures;
    }

    std::vector<LoadedCapture> selected;
    selected.reserve(static_cast<std::size_t>(requested_count));
    const auto denominator = requested_count - 1;
    for (int index = 0; index < requested_count; ++index)
    {
        const auto source_index = denominator == 0
            ? 0
            : static_cast<int>(RoundToEven(
                static_cast<double>(index) * static_cast<double>(captures.size() - 1U) /
                static_cast<double>(denominator)));
        selected.push_back(captures[static_cast<std::size_t>(source_index)]);
    }

    return selected;
}

void WriteDepthValidationSummary(const ValidationSummary& summary, const fs::path& output_directory)
{
    std::ofstream stream(output_directory / "validation_summary.txt", std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to write validation_summary.txt.");
    }

    stream << "Validation summary\n";
    stream << "label=" << summary.label << "\n";
    stream << "images_checked=" << summary.images_checked << "\n";
    stream << "points_checked=" << summary.points_checked << "\n";
    stream << "positive_depth_count=" << summary.positive_depth_count << "\n";
    stream << "in_frame_count=" << summary.in_frame_count << "\n";
    stream << "depth_compared_count=" << summary.depth_compared_count << "\n";
    stream << "mean_depth_abs_error=" << ToRoundTrip(summary.mean_depth_abs_error) << "\n";
    stream << "max_depth_abs_error=" << ToRoundTrip(summary.max_depth_abs_error) << "\n";
}

ValidationSummary RunDepthValidation(
    const LoadedDataset& dataset,
    const std::vector<Vector3d>& world_points,
    const std::string& label,
    const fs::path& output_directory,
    const int image_count,
    const int point_stride)
{
    const auto captures = SelectValidationCaptures(dataset.captures, image_count);
    auto points_checked = 0;
    auto positive_depth_count = 0;
    auto in_frame_count = 0;
    auto depth_compared_count = 0;
    auto depth_abs_error_sum = 0.0;
    auto max_depth_abs_error = 0.0;

    for (const auto& capture : captures)
    {
        if (!capture.depth_info.has_value())
        {
            continue;
        }

        const auto depth_values = LoadDepthMap(*capture.depth_info);
        for (std::size_t point_index = 0; point_index < world_points.size(); point_index += static_cast<std::size_t>(point_stride))
        {
            ++points_checked;
            const auto projection = ProjectWorldPoint(capture, world_points[point_index]);
            if (!projection.has_value())
            {
                continue;
            }

            ++positive_depth_count;
            if (projection->pixel_x < 0.0 ||
                projection->pixel_x >= static_cast<double>(capture.width) ||
                projection->pixel_y < 0.0 ||
                projection->pixel_y >= static_cast<double>(capture.height))
            {
                continue;
            }

            ++in_frame_count;
            const auto depth_pixel = MapImagePixelToDepthPixel(capture, projection->pixel_x, projection->pixel_y);
            if (!depth_pixel.has_value())
            {
                continue;
            }

            const auto depth_meters = depth_values[
                static_cast<std::size_t>(depth_pixel->y) * static_cast<std::size_t>(capture.depth_info->width) +
                static_cast<std::size_t>(depth_pixel->x)];
            if (!(depth_meters > 0.0F) || !std::isfinite(depth_meters))
            {
                continue;
            }

            ++depth_compared_count;
            const auto abs_error = std::abs(static_cast<double>(depth_meters) - projection->depth_meters);
            depth_abs_error_sum += abs_error;
            max_depth_abs_error = std::max(max_depth_abs_error, abs_error);
        }
    }

    const auto summary = ValidationSummary{
        label,
        static_cast<int>(captures.size()),
        points_checked,
        positive_depth_count,
        in_frame_count,
        depth_compared_count,
        depth_compared_count == 0 ? 0.0 : (depth_abs_error_sum / static_cast<double>(depth_compared_count)),
        max_depth_abs_error};
    WriteDepthValidationSummary(summary, output_directory);
    return summary;
}

void WriteReprojectionSummary(const ReprojectionSummary& summary, const fs::path& output_directory)
{
    std::ofstream stream(output_directory / "reprojection_summary.txt", std::ios::binary);
    if (!stream)
    {
        throw std::runtime_error("Failed to write reprojection_summary.txt.");
    }

    stream << "Reprojection summary\n";
    stream << "points_checked=" << summary.point_count << "\n";
    stream << "observations_checked=" << summary.observation_count << "\n";
    stream << "mean_pixel_error=" << ToRoundTrip(summary.mean_pixel_error) << "\n";
    stream << "max_pixel_error=" << ToRoundTrip(summary.max_pixel_error) << "\n";
}

ReprojectionSummary RunReprojectionValidation(
    const LoadedDataset& dataset,
    const SparseReconstruction& sparse_reconstruction,
    const fs::path& output_directory)
{
    std::map<int, const LoadedCapture*> captures_by_id;
    std::map<int, const SparseImagePoints*> observations_by_image;
    for (const auto& capture : dataset.captures)
    {
        captures_by_id[capture.index] = &capture;
    }
    for (const auto& image : sparse_reconstruction.image_points)
    {
        observations_by_image[image.image_id] = &image;
    }

    auto total_pixel_error = 0.0;
    auto max_pixel_error = 0.0;
    auto observation_count = 0;
    for (const auto& point : sparse_reconstruction.points3d)
    {
        for (const auto& track_element : point.track)
        {
            const auto capture_it = captures_by_id.find(track_element.image_id);
            const auto observation_it = observations_by_image.find(track_element.image_id);
            if (capture_it == captures_by_id.end() ||
                observation_it == observations_by_image.end() ||
                track_element.point2d_index < 0 ||
                static_cast<std::size_t>(track_element.point2d_index) >= observation_it->second->observations.size())
            {
                continue;
            }

            const auto& capture = *capture_it->second;
            const auto& observation = observation_it->second->observations[static_cast<std::size_t>(track_element.point2d_index)];
            const auto projection = ProjectWorldPoint(capture, point.position);
            if (!projection.has_value())
            {
                continue;
            }

            const auto dx = projection->pixel_x - observation.x;
            const auto dy = projection->pixel_y - observation.y;
            const auto pixel_error = std::sqrt((dx * dx) + (dy * dy));
            total_pixel_error += pixel_error;
            max_pixel_error = std::max(max_pixel_error, pixel_error);
            ++observation_count;
        }
    }

    const auto summary = ReprojectionSummary{
        static_cast<int>(sparse_reconstruction.points3d.size()),
        observation_count,
        observation_count == 0 ? 0.0 : (total_pixel_error / static_cast<double>(observation_count)),
        max_pixel_error};
    WriteReprojectionSummary(summary, output_directory);
    return summary;
}

}  // namespace

bool ColmapModelValidationSummary::IsMatch() const
{
    return text_counts.cameras == binary_counts.cameras &&
        text_counts.images == binary_counts.images &&
        text_counts.points3d == binary_counts.points3d &&
        text_counts.image_observations == binary_counts.image_observations &&
        text_counts.point_track_elements == binary_counts.point_track_elements;
}

HelpRequested::HelpRequested()
    : std::runtime_error("Help requested.")
{
}

CliOptions CreateDefaultOptions(
    const fs::path& input_directory,
    const fs::path& output_directory,
    const bool validate_exported_model,
    const bool copy_images)
{
    return {
        fs::absolute(input_directory),
        fs::absolute(output_directory),
        std::nullopt,
        false,
        copy_images,
        validate_exported_model,
        5,
        256,
        {4, 6, 4, 0.12, 2.0, 0.50, 0.20}};
}

CliOptions ParseCliOptions(const std::vector<std::string>& args)
{
    fs::path input_directory;
    fs::path output_directory;
    std::optional<fs::path> align_to_colmap_model_path;
    auto snap_shared_camera_poses_to_reference = false;
    auto copy_images = true;
    auto validate_exported_model = false;
    auto validation_image_count = 5;
    auto validation_point_stride = 256;
    auto depth_bucket_size = 4;
    auto max_observations_per_track = 6;
    auto min_capture_gap = 4;
    auto min_baseline_meters = 0.12;
    auto min_angular_separation_degrees = 2.0;
    auto base_depth_tolerance_meters = 0.50;
    auto relative_depth_tolerance = 0.20;

    for (std::size_t index = 0; index < args.size(); ++index)
    {
        const auto& arg = args[index];
        if (arg == "--input" || arg == "-i")
        {
            input_directory = RequireValue(args, index, arg);
        }
        else if (arg == "--output" || arg == "-o")
        {
            output_directory = RequireValue(args, index, arg);
        }
        else if (arg == "--align-to-colmap-model")
        {
            align_to_colmap_model_path = RequireValue(args, index, arg);
        }
        else if (arg == "--snap-shared-cameras-to-reference")
        {
            snap_shared_camera_poses_to_reference = true;
        }
        else if (arg == "--no-copy-images")
        {
            copy_images = false;
        }
        else if (arg == "--validate-exported-model")
        {
            validate_exported_model = true;
        }
        else if (arg == "--validation-images")
        {
            validation_image_count = ParsePositiveInt(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--validation-stride")
        {
            validation_point_stride = ParsePositiveInt(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--depth-bucket-size")
        {
            depth_bucket_size = ParsePositiveInt(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--max-observations-per-track")
        {
            max_observations_per_track = ParsePositiveInt(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--min-capture-gap")
        {
            min_capture_gap = ParseNonNegativeInt(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--min-baseline-m")
        {
            min_baseline_meters = ParsePositiveDouble(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--min-angle-deg")
        {
            min_angular_separation_degrees = ParsePositiveDouble(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--base-depth-tolerance-m")
        {
            base_depth_tolerance_meters = ParsePositiveDouble(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--relative-depth-tolerance")
        {
            relative_depth_tolerance = ParsePositiveDouble(RequireValue(args, index, arg), arg);
        }
        else if (arg == "--help" || arg == "-h")
        {
            throw HelpRequested();
        }
        else
        {
            throw std::invalid_argument("Unknown argument '" + arg + "'.");
        }
    }

    if (input_directory.empty())
    {
        throw std::invalid_argument("Missing required argument --input.");
    }

    if (output_directory.empty())
    {
        output_directory = fs::current_path() / "output";
    }

    return {
        fs::absolute(input_directory),
        fs::absolute(output_directory),
        align_to_colmap_model_path.has_value() ? std::optional<fs::path>(fs::absolute(*align_to_colmap_model_path)) : std::nullopt,
        snap_shared_camera_poses_to_reference,
        copy_images,
        validate_exported_model,
        validation_image_count,
        validation_point_stride,
        {depth_bucket_size,
         max_observations_per_track,
         min_capture_gap,
         min_baseline_meters,
         min_angular_separation_degrees,
         base_depth_tolerance_meters,
         relative_depth_tolerance}};
}

std::string GetUsageText()
{
    return
        "Usage:\n"
        "  sk2colmap --input <SplatkingFolder> [--output <OutputFolder>] "
        "[--align-to-colmap-model <ColmapModelPath>] "
        "[--snap-shared-cameras-to-reference] [--no-copy-images] [--validate-exported-model] "
        "[--validation-images <N>] [--validation-stride <N>] [--depth-bucket-size <N>] "
        "[--max-observations-per-track <N>] [--min-capture-gap <N>] [--min-baseline-m <meters>] "
        "[--min-angle-deg <degrees>] [--base-depth-tolerance-m <meters>] "
        "[--relative-depth-tolerance <ratio>]\n\n"
        "Example:\n"
        "  sk2colmap --input ./LidarSeries --output ./output";
}

ConversionResult RunConversion(
    const CliOptions& options,
    const std::function<void(const std::string&)>& log)
{
    const auto log_message = [&](const std::string& message) {
        if (log)
        {
            log(message);
        }
    };

    log_message("Loading dataset from " + options.input_directory.string());
    auto dataset = LoadDataset(options.input_directory);
    log_message(
        "Loaded " +
        std::to_string(dataset.captures.size()) +
        " captures and " +
        std::to_string(dataset.world_points.size()) +
        " LiDAR points.");

    auto reconstruction = BuildSparseReconstruction(dataset, options.refinement);
    auto total_observations = 0;
    for (const auto& image : reconstruction.sparse_reconstruction.image_points)
    {
        total_observations += static_cast<int>(image.observations.size());
    }

    auto average_track_length = reconstruction.sparse_reconstruction.points3d.empty()
        ? 0.0
        : static_cast<double>(total_observations) /
            static_cast<double>(reconstruction.sparse_reconstruction.points3d.size());
    log_message(
        "Synthesized " +
        std::to_string(reconstruction.sparse_reconstruction.points3d.size()) +
        " sparse LiDAR points from " +
        std::to_string(reconstruction.dense_aligned_points.size()) +
        " raw dense samples with " +
        std::to_string(total_observations) +
        " 2D observations (avg track length " +
        ToFixed(average_track_length, 2) +
        ").");

    auto export_dataset = dataset;
    auto export_reconstruction = reconstruction;
    std::optional<AlignmentSummary> alignment_summary;
    if (options.align_to_colmap_model_path.has_value())
    {
        auto alignment = ApplyReferenceAlignment(
            dataset,
            reconstruction,
            *options.align_to_colmap_model_path,
            options.snap_shared_camera_poses_to_reference);
        export_dataset = std::move(alignment.dataset);
        export_reconstruction = std::move(alignment.reconstruction);
        alignment_summary = alignment.summary;

        total_observations = 0;
        for (const auto& image : export_reconstruction.sparse_reconstruction.image_points)
        {
            total_observations += static_cast<int>(image.observations.size());
        }

        average_track_length = export_reconstruction.sparse_reconstruction.points3d.empty()
            ? 0.0
            : static_cast<double>(total_observations) /
                static_cast<double>(export_reconstruction.sparse_reconstruction.points3d.size());

        log_message(
            "Aligned to reference COLMAP model with " +
            std::to_string(alignment_summary->matched_image_count) +
            " shared images, scale=" +
            ToFixed(alignment_summary->scale, 6) +
            ", mean center residual=" +
            ToFixed(alignment_summary->mean_center_residual_meters, 4) +
            " m, max=" +
            ToFixed(alignment_summary->max_center_residual_meters, 4) +
            " m.");
        if (alignment_summary->shared_camera_poses_snapped_to_reference)
        {
            log_message(
                "Snapped " +
                std::to_string(alignment_summary->snapped_image_count) +
                " shared camera poses exactly to the reference model; reprojection consistency with LiDAR may degrade.");
        }
        log_message(
            "Transformed sparse model into reference frame with " +
            std::to_string(export_reconstruction.sparse_reconstruction.points3d.size()) +
            " LiDAR points, " +
            std::to_string(total_observations) +
            " 2D observations and avg track length " +
            ToFixed(average_track_length, 2) +
            ".");
    }

    WriteColmapOutput(export_dataset, export_reconstruction, options.output_directory, options.copy_images);
    if (alignment_summary.has_value())
    {
        WriteAlignmentSummary(*alignment_summary, options.output_directory);
    }
    log_message("COLMAP export written to " + options.output_directory.string());

    std::optional<ColmapModelValidationSummary> model_validation;
    if (options.validate_exported_model)
    {
        model_validation = ValidateExportedModel(options.output_directory);
        log_message(
            "Model validation: match=" +
            BoolText(model_validation->IsMatch()) +
            ", cameras=" +
            std::to_string(model_validation->binary_counts.cameras) +
            ", images=" +
            std::to_string(model_validation->binary_counts.images) +
            ", points3D=" +
            std::to_string(model_validation->binary_counts.points3d) +
            ".");
    }

    const auto camera_pose_validation = RunCameraPoseValidation(export_dataset, options.output_directory);
    const auto validation = RunDepthValidation(
        dataset,
        dataset.world_points,
        "raw_session_lidar_cloud",
        options.output_directory,
        options.validation_image_count,
        options.validation_point_stride);
    const auto reprojection = RunReprojectionValidation(
        export_dataset,
        export_reconstruction.sparse_reconstruction,
        options.output_directory);

    log_message(
        "Camera pose validation: " +
        std::to_string(camera_pose_validation.camera_count) +
        " cameras checked, mean center error = " +
        ToRoundTrip(camera_pose_validation.mean_center_error_meters) +
        " m, max = " +
        ToRoundTrip(camera_pose_validation.max_center_error_meters) +
        " m.");
    log_message(
        "Validation (" +
        validation.label +
        "): " +
        std::to_string(validation.depth_compared_count) +
        " depth samples compared, mean |dz| = " +
        ToFixed(validation.mean_depth_abs_error, 4) +
        " m, max |dz| = " +
        ToFixed(validation.max_depth_abs_error, 4) +
        " m.");
    log_message(
        "Sparse reprojection: " +
        std::to_string(reprojection.observation_count) +
        " observations compared, mean pixel error = " +
        ToFixed(reprojection.mean_pixel_error, 2) +
        "px, max = " +
        ToFixed(reprojection.max_pixel_error, 2) +
        "px.");

    return {
        options.output_directory,
        static_cast<int>(dataset.captures.size()),
        static_cast<int>(export_reconstruction.dense_aligned_points.size()),
        static_cast<int>(export_reconstruction.sparse_reconstruction.points3d.size()),
        total_observations,
        average_track_length,
        alignment_summary,
        validation,
        model_validation};
}

}  // namespace sk2cm
