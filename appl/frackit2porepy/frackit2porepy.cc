#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <variant>
#include <array>
#include <cmath>
#include <type_traits>
#include <memory>

#include <frackit/common/id.hh>
#include <frackit/io/gmshwriter.hh>

#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/sampling/disksampler.hh>
#include <frackit/sampling/quadrilateralsampler.hh>
#include <frackit/sampling/multigeometrysampler.hh>
#include <frackit/entitynetwork/constraintsmatrix.hh>
#include <frackit/entitynetwork/multigeometryentityset.hh>
#include <frackit/sampling/status.hh>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/box.hh>

#include <frackit/entitynetwork/constraints.hh>
#include <frackit/entitynetwork/networkbuilder.hh>

#include "toml.hpp" // Third-party TOML library

using DiskSampler_t = Frackit::DiskSampler<double>;
using QuadrilateralSampler_t = Frackit::QuadrilateralSampler<3, double>;

namespace
{

    struct NormalParams
    {
        double mean = 0.0;
        double stddev = 0.0;
    };

    struct UniformParams
    {
        double min = 0.0;
        double max = 1.0;
    };

    struct DomainConfig
    {
        double xmin = 0, ymin = 0, zmin = 0, xmax = 100, ymax = 100, zmax = 100;
    };

    struct DiskSamplerConfig
    {
        NormalParams major_axis_length{50.0, 10.0};
        NormalParams minor_axis_length{50.0, 5.0};
        NormalParams rot_x_deg{0.0, 7.5};
        NormalParams rot_y_deg{0.0, 7.5};
        NormalParams rot_z_deg{0.0, 7.5};
    };

    struct QuadSamplerConfig
    {
        NormalParams strike_deg{0.0, 7.5};
        NormalParams dip_deg{0.0, 7.5};
        UniformParams strike_length{50.0, 10.0};
        UniformParams dip_length{50.0, 5.0};
    };

    struct ConstraintsConfig
    {
        double min_distance = 0.05;
        double min_intersecting_angle_deg_self = 30.0;
        double min_intersecting_angle_deg_other = 40.0;
        double min_intersection_magnitude = 0.05;
        double min_intersection_distance = 0.05;
    };

    using AnySamplerConfig = std::variant<DiskSamplerConfig, QuadSamplerConfig>;

    struct FamilyConfig
    {
        int target_num = 5;
        std::string type = "disks";                     // "disks" or "quads"
        AnySamplerConfig sampler = DiskSamplerConfig{}; // default matches default type
    };

    struct OutputConfig
    {
        std::string disks_csv = "disks.csv";
        std::string quads_csv = "quads.csv";
        std::string families_csv = "families.csv";
    };

    struct Config
    {
        DomainConfig domain{};
        DomainConfig subdomain{};
        ConstraintsConfig constraints{};
        std::vector<FamilyConfig> families{};
        OutputConfig output{};
    };

    static void printConfig(const Config &cfg)
    {
        std::cout << "\n[config] domain: [" << cfg.domain.xmin << ", " << cfg.domain.ymin << ", " << cfg.domain.zmin
                  << "] -> [" << cfg.domain.xmax << ", " << cfg.domain.ymax << ", " << cfg.domain.zmax << "]\n";
        std::cout << "\n[config] subdomain: [" << cfg.subdomain.xmin << ", " << cfg.subdomain.ymin << ", " << cfg.subdomain.zmin
                  << "] -> [" << cfg.subdomain.xmax << ", " << cfg.subdomain.ymax << ", " << cfg.subdomain.zmax << "]\n";
        if (cfg.families.empty())
            std::cout << "[config] sampler.num=0\n";
        else
        {
            std::cout << "[config] sampler.num=" << cfg.families.size() << "\n";
            for (std::size_t i = 0; i < cfg.families.size(); ++i)
            {
                const auto idx = i + 1;
                const auto &f = cfg.families[i];
                std::cout << "[config] sampler." << idx << ".type=\"" << f.type << "\"\n";
                std::cout << "[config] sampler." << idx << ".target_num=" << f.target_num << "\n";
                std::visit(
                    [&](const auto &s)
                    {
                        using T = std::decay_t<decltype(s)>;
                        if constexpr (std::is_same_v<T, DiskSamplerConfig>)
                        {
                            std::cout << "[config] sampler." << idx << ".major_axis_length.mean=" << s.major_axis_length.mean << "\n";
                            std::cout << "[config] sampler." << idx << ".major_axis_length.stddev=" << s.major_axis_length.stddev << "\n";
                            std::cout << "[config] sampler." << idx << ".minor_axis_length.mean=" << s.minor_axis_length.mean << "\n";
                            std::cout << "[config] sampler." << idx << ".minor_axis_length.stddev=" << s.minor_axis_length.stddev << "\n";
                            std::cout << "[config] sampler." << idx << ".rot_x_deg.mean=" << s.rot_x_deg.mean << "\n";
                            std::cout << "[config] sampler." << idx << ".rot_x_deg.stddev=" << s.rot_x_deg.stddev << "\n";
                            std::cout << "[config] sampler." << idx << ".rot_y_deg.mean=" << s.rot_y_deg.mean << "\n";
                            std::cout << "[config] sampler." << idx << ".rot_y_deg.stddev=" << s.rot_y_deg.stddev << "\n";
                            std::cout << "[config] sampler." << idx << ".rot_z_deg.mean=" << s.rot_z_deg.mean << "\n";
                            std::cout << "[config] sampler." << idx << ".rot_z_deg.stddev=" << s.rot_z_deg.stddev << "\n\n";
                        }
                        else if constexpr (std::is_same_v<T, QuadSamplerConfig>)
                        {
                            std::cout << "[config] sampler." << idx << ".strike_deg.mean=" << s.strike_deg.mean << "\n";
                            std::cout << "[config] sampler." << idx << ".strike_deg.stddev=" << s.strike_deg.stddev << "\n";
                            std::cout << "[config] sampler." << idx << ".dip_deg.mean=" << s.dip_deg.mean << "\n";
                            std::cout << "[config] sampler." << idx << ".dip_deg.stddev=" << s.dip_deg.stddev << "\n";
                            std::cout << "[config] sampler." << idx << ".strike_length.min=" << s.strike_length.min << "\n";
                            std::cout << "[config] sampler." << idx << ".strike_length.max=" << s.strike_length.max << "\n";
                            std::cout << "[config] sampler." << idx << ".dip_length.min=" << s.dip_length.min << "\n";
                            std::cout << "[config] sampler." << idx << ".dip_length.max=" << s.dip_length.max << "\n\n";
                        }
                    },
                    f.sampler);
            }
        }

        std::cout << "[config] constraints.min_distance=" << cfg.constraints.min_distance << "\n";
        std::cout << "[config] constraints.min_intersecting_angle_deg_self=" << cfg.constraints.min_intersecting_angle_deg_self << "\n";
        std::cout << "[config] constraints.min_intersecting_angle_deg_other=" << cfg.constraints.min_intersecting_angle_deg_other << "\n";
        std::cout << "[config] constraints.min_intersection_magnitude=" << cfg.constraints.min_intersection_magnitude << "\n";
        std::cout << "[config] constraints.min_intersection_distance=" << cfg.constraints.min_intersection_distance << "\n";
        std::cout << "[config] output.disks_csv=\"" << cfg.output.disks_csv << "\"\n";
        std::cout << "[config] output.families_csv=\"" << cfg.output.families_csv << "\"\n";
    }

    static void tryLoadToml(const std::string &path, Config &cfg)
    {
        if (!std::ifstream(path).good())
        {
            std::cout << "[config] Note: config file not found (using defaults): " << path << "\n";
            return;
        }

        try
        {
            const auto tbl = toml::parse_file(path);

            if (const auto *domain = tbl["domain"].as_table())
            {
                cfg.domain.xmin = (*domain)["xmin"].value_or(cfg.domain.xmin);
                cfg.domain.ymin = (*domain)["ymin"].value_or(cfg.domain.ymin);
                cfg.domain.zmin = (*domain)["zmin"].value_or(cfg.domain.zmin);
                cfg.domain.xmax = (*domain)["xmax"].value_or(cfg.domain.xmax);
                cfg.domain.ymax = (*domain)["ymax"].value_or(cfg.domain.ymax);
                cfg.domain.zmax = (*domain)["zmax"].value_or(cfg.domain.zmax);
            }
            if (const auto *subdomain = tbl["subdomain"].as_table())
            {
                cfg.subdomain.xmin = (*subdomain)["xmin"].value_or(cfg.subdomain.xmin);
                cfg.subdomain.ymin = (*subdomain)["ymin"].value_or(cfg.subdomain.ymin);
                cfg.subdomain.zmin = (*subdomain)["zmin"].value_or(cfg.subdomain.zmin);
                cfg.subdomain.xmax = (*subdomain)["xmax"].value_or(cfg.subdomain.xmax);
                cfg.subdomain.ymax = (*subdomain)["ymax"].value_or(cfg.subdomain.ymax);
                cfg.subdomain.zmax = (*subdomain)["zmax"].value_or(cfg.subdomain.zmax);
            }

            if (!(cfg.domain.xmin <= cfg.subdomain.xmin && cfg.domain.ymin <= cfg.subdomain.ymin && cfg.domain.zmin <= cfg.subdomain.zmin &&
                  cfg.domain.xmax >= cfg.subdomain.xmax && cfg.domain.ymax >= cfg.subdomain.ymax && cfg.domain.zmax >= cfg.subdomain.zmax))
            {
                throw std::runtime_error("Subdomain must be contained within the domain.");
            }

            auto loadDiskSamplerFromTable = [](const toml::table &st, DiskSamplerConfig &s)
            {
                if (const auto *maj = st["major_axis_length"].as_table())
                {
                    s.major_axis_length.mean = (*maj)["mean"].value_or(s.major_axis_length.mean);
                    s.major_axis_length.stddev = (*maj)["stddev"].value_or(s.major_axis_length.stddev);
                }
                if (const auto *min = st["minor_axis_length"].as_table())
                {
                    s.minor_axis_length.mean = (*min)["mean"].value_or(s.minor_axis_length.mean);
                    s.minor_axis_length.stddev = (*min)["stddev"].value_or(s.minor_axis_length.stddev);
                }

                if (const auto *rot = st["rotation_deg"].as_table())
                {
                    if (const auto *rx = (*rot)["x"].as_table())
                    {
                        s.rot_x_deg.mean = (*rx)["mean"].value_or(s.rot_x_deg.mean);
                        s.rot_x_deg.stddev = (*rx)["stddev"].value_or(s.rot_x_deg.stddev);
                    }
                    if (const auto *ry = (*rot)["y"].as_table())
                    {
                        s.rot_y_deg.mean = (*ry)["mean"].value_or(s.rot_y_deg.mean);
                        s.rot_y_deg.stddev = (*ry)["stddev"].value_or(s.rot_y_deg.stddev);
                    }
                    if (const auto *rz = (*rot)["z"].as_table())
                    {
                        s.rot_z_deg.mean = (*rz)["mean"].value_or(s.rot_z_deg.mean);
                        s.rot_z_deg.stddev = (*rz)["stddev"].value_or(s.rot_z_deg.stddev);
                    }
                }
            };

            auto loadQuadSamplerFromTable = [](const toml::table &st, QuadSamplerConfig &s)
            {
                if (const auto *strike = st["strike_deg"].as_table())
                {
                    s.strike_deg.mean = (*strike)["mean"].value_or(s.strike_deg.mean);
                    s.strike_deg.stddev = (*strike)["stddev"].value_or(s.strike_deg.stddev);
                }
                if (const auto *dip = st["dip_deg"].as_table())
                {
                    s.dip_deg.mean = (*dip)["mean"].value_or(s.dip_deg.mean);
                    s.dip_deg.stddev = (*dip)["stddev"].value_or(s.dip_deg.stddev);
                }
                if (const auto *sl = st["strike_length"].as_table())
                {
                    s.strike_length.min = (*sl)["min"].value_or(s.strike_length.min);
                    s.strike_length.max = (*sl)["max"].value_or(s.strike_length.max);
                }
                if (const auto *dl = st["dip_length"].as_table())
                {
                    s.dip_length.min = (*dl)["min"].value_or(s.dip_length.min);
                    s.dip_length.max = (*dl)["max"].value_or(s.dip_length.max);
                }
            };

            if (const auto *sampler = tbl["sampler"].as_table())
            {
                const int num = (*sampler)["num"].value_or(0);
                if (num > 0)
                {
                    cfg.families.clear();
                    cfg.families.resize(static_cast<std::size_t>(num));

                    for (int i = 1; i <= num; ++i)
                    {
                        auto &f = cfg.families[static_cast<std::size_t>(i - 1)];
                        const auto key = std::to_string(i);

                        const auto *ft = (*sampler)[key].as_table();
                        if (!ft)
                            continue;

                        f.target_num = (*ft)["target_num"].value_or(f.target_num);
                        f.type = (*ft)["type"].value_or(f.type);

                        if (f.type == "disks")
                        {
                            DiskSamplerConfig s{};
                            loadDiskSamplerFromTable(*ft, s);
                            f.sampler = s;
                        }
                        else if (f.type == "quads")
                        {
                            QuadSamplerConfig s{};
                            loadQuadSamplerFromTable(*ft, s);
                            f.sampler = s;
                        }
                        else
                        {
                            throw std::runtime_error("Unsupported sampler type in TOML: " + f.type);
                        }
                    }
                }
            }

            if (const auto *c = tbl["constraints"].as_table())
            {
                cfg.constraints.min_distance = (*c)["min_distance"].value_or(cfg.constraints.min_distance);
                cfg.constraints.min_intersecting_angle_deg_self =
                    (*c)["min_intersecting_angle_deg_self"].value_or(cfg.constraints.min_intersecting_angle_deg_self);
                cfg.constraints.min_intersecting_angle_deg_other =
                    (*c)["min_intersecting_angle_deg_other"].value_or(cfg.constraints.min_intersecting_angle_deg_other);
                cfg.constraints.min_intersection_magnitude =
                    (*c)["min_intersection_magnitude"].value_or(cfg.constraints.min_intersection_magnitude);
                cfg.constraints.min_intersection_distance =
                    (*c)["min_intersection_distance"].value_or(cfg.constraints.min_intersection_distance);
            }

            if (const auto *o = tbl["output"].as_table())
            {
                cfg.output.disks_csv = (*o)["disks_csv"].value_or(cfg.output.disks_csv);
                cfg.output.quads_csv = (*o)["quads_csv"].value_or(cfg.output.quads_csv);
                cfg.output.families_csv = (*o)["families_csv"].value_or(cfg.output.families_csv);
            }
        }
        catch (const toml::parse_error &err)
        {
            std::cerr << "[config] TOML parse error in " << path << ": " << err.description() << "\n";
        }
    }

} // namespace

// Helpers for conversion to PorePy format

using Vec3 = std::array<double, 3>;

auto dot = [](const Vec3 &a, const Vec3 &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};

auto cross = [](const Vec3 &a, const Vec3 &b) -> Vec3
{
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]};
};

auto norm = [&](const Vec3 &v)
{
    return std::sqrt(dot(v, v));
};

auto normalize = [&](const Vec3 &v) -> Vec3
{
    const double n = norm(v);
    if (n < 1e-14)
        throw std::runtime_error("normalize(): zero-length vector");
    return {v[0] / n, v[1] / n, v[2] / n};
};

auto wrap_0_2pi = [](double a)
{
    constexpr double two_pi = 2.0 * M_PI;
    a = std::fmod(a, two_pi);
    if (a < 0.0)
        a += two_pi;
    return a;
};

// Rodrigues rotation of a vector around an axis
auto rotate_about_axis = [&](const Vec3 &v, const Vec3 &axis_raw, double theta) -> Vec3
{
    const Vec3 u = normalize(axis_raw);
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const Vec3 uxv = cross(u, v);
    const double udv = dot(u, v);

    return {
        v[0] * c + uxv[0] * s + u[0] * udv * (1.0 - c),
        v[1] * c + uxv[1] * s + u[1] * udv * (1.0 - c),
        v[2] * c + uxv[2] * s + u[2] * udv * (1.0 - c)};
};

auto toPorePyAngles = [&](Vec3 major_axis_dir, Vec3 normal_dir)
{
    // Normalize inputs
    Vec3 a = normalize(major_axis_dir);
    Vec3 n = normalize(normal_dir);

    // Enforce upward normal for consistency
    if (n[2] < 0.0)
    {
        n = {-n[0], -n[1], -n[2]};
        a = {-a[0], -a[1], -a[2]}; // keep in-plane orientation consistent
    }

    // Dip: rotation around strike axis (radians)
    const double dip = std::acos(std::max(-1.0, std::min(1.0, n[2])));

    // Strike direction = cross(k, n) (horizontal), k=(0,0,1)
    Vec3 strike_dir = {-n[1], n[0], 0.0};
    const double s_norm = norm(strike_dir);

    double strike = 0.0;
    if (s_norm > 1e-12)
    {
        strike_dir = {strike_dir[0] / s_norm, strike_dir[1] / s_norm, 0.0};

        // Canonicalize to reduce flips (optional)
        if (strike_dir[1] < 0.0)
        {
            strike_dir = {-strike_dir[0], -strike_dir[1], 0.0};
        }

        // PorePy strike angle: clockwise from +y (north)
        // strike_dir = (sin(strike), cos(strike), 0)
        strike = std::atan2(strike_dir[0], strike_dir[1]);
        strike = wrap_0_2pi(strike);
    }

    // Strike axis as PorePy uses it
    const Vec3 strike_axis = {std::sin(strike), std::cos(strike), 0.0};

    // Undo the dip tilt to recover the pre-strike-dip in-plane major axis direction
    Vec3 a0 = a;
    if (dip > 1e-14)
    {
        a0 = rotate_about_axis(a, strike_axis, -dip);
    }

    // Major axis angle: rotation from +x in the XY plane before strike-dip
    const double major_axis_angle = wrap_0_2pi(std::atan2(a0[1], a0[0]));

    return std::tuple<double, double, double>{major_axis_angle, strike, dip};
};

//! create a network of 3d quadrilaterals
int main(int argc, char **argv)
{
    using namespace Frackit;

    // Create a configuration object
    Config cfg;
    tryLoadToml("/frackit/shared/config.toml", cfg);
    printConfig(cfg);

    // We consider 3d space here
    static constexpr int worldDimension = 3;

    // Define the type used for coordinates
    using ctype = double;

    // Define a domain (here: unit cube) in which the subdomain lives.
    // Boxes are created by providing xmin, ymin, zmin and xmax, ymax and zmax in constructor.
    Box<ctype> domain(cfg.domain.xmin, cfg.domain.ymin, cfg.domain.zmin,
                      cfg.domain.xmax, cfg.domain.ymax, cfg.domain.zmax);

    // Define a subdomain (here: unit cube) in which the fractures should be created.
    // Boxes are created by providing xmin, ymin, zmin and xmax, ymax and zmax in constructor.
    Box<ctype> subdomain(cfg.subdomain.xmin, cfg.subdomain.ymin, cfg.subdomain.zmin,
                         cfg.subdomain.xmax, cfg.subdomain.ymax, cfg.subdomain.zmax);

    // We now create a sampler instance that uniformly samples points within this box.
    // These points will be used as the quadrilateral centers.
    auto pointSampler = makeUniformPointSampler(subdomain);

    // Sampler class for quadrilaterals. Per default, this uses uniform distributions
    // for all parameters defining the quadrilaterals. Quadrilateral samplers require
    // distributions for strike angle, dip angle, edge length (see class description
    // for more details).
    using NormalDistro = std::normal_distribution<ctype>;
    using UniformDistro = std::uniform_real_distribution<ctype>;

    auto makeDiskSamplerFromCfg = [&](const DiskSamplerConfig &sc)
    {
        return DiskSampler_t(
            pointSampler,
            NormalDistro(sc.major_axis_length.mean, sc.major_axis_length.stddev),
            NormalDistro(sc.minor_axis_length.mean, sc.minor_axis_length.stddev),
            NormalDistro(toRadians(sc.rot_x_deg.mean), toRadians(sc.rot_x_deg.stddev)),
            NormalDistro(toRadians(sc.rot_y_deg.mean), toRadians(sc.rot_y_deg.stddev)),
            NormalDistro(toRadians(sc.rot_z_deg.mean), toRadians(sc.rot_z_deg.stddev)));
    };

    auto makeQuadSamplerFromCfg = [&](const QuadSamplerConfig &sc)
    {
        return QuadrilateralSampler_t(
            pointSampler,
            NormalDistro(toRadians(sc.strike_deg.mean), toRadians(sc.strike_deg.stddev)),
            NormalDistro(toRadians(sc.dip_deg.mean), toRadians(sc.dip_deg.stddev)),
            UniformDistro(sc.strike_length.min, sc.strike_length.max),
            UniformDistro(sc.dip_length.min, sc.dip_length.max));
    };

    // Central sampler
    using Disk = Disk<ctype>;
    using Quad = Quadrilateral<ctype, worldDimension>;
    MultiGeometrySampler<Disk, Quad> multiSampler;

    // Central status
    SamplingStatus status;
    for (std::size_t i = 0; i < cfg.families.size(); ++i)
    {
        const auto &fam = cfg.families[i];
        const Id familyId(static_cast<int>(i + 1));

        status.setTargetCount(familyId, fam.target_num);

        if (fam.type == "disks")
        {
            const auto sc = std::get<DiskSamplerConfig>(fam.sampler);
            auto diskSampler = makeDiskSamplerFromCfg(sc);
            multiSampler.addGeometrySampler(diskSampler, familyId);
        }
        else if (fam.type == "quads")
        {
            const auto sc = std::get<QuadSamplerConfig>(fam.sampler);
            auto quadSampler = makeQuadSamplerFromCfg(sc);
            multiSampler.addGeometrySampler(quadSampler, familyId);
        }
        else
            throw std::runtime_error("Unsupported family type: " + fam.type);
    }

    // We want to enforce some constraints.
    // In particular, for entities of the same set we want a minimum spacing
    // distance, and the entities must not intersect in angles
    // less than given value. Moreover, if they intersect, we don't want intersection
    // edges whose length is smaller than given value, and, the intersection should not
    // be too close to the boundary of one of two intersecting entities.
    using Constraints = EntityNetworkConstraints<ctype>;
    Constraints cSelf;
    cSelf.setMinDistance(cfg.constraints.min_distance);
    cSelf.setMinIntersectingAngle(toRadians(cfg.constraints.min_intersecting_angle_deg_self));
    cSelf.setMinIntersectionMagnitude(cfg.constraints.min_intersection_magnitude);
    cSelf.setMinIntersectionDistance(cfg.constraints.min_intersection_distance);

    Constraints cOther = cSelf;
    cOther.setMinIntersectingAngle(toRadians(cfg.constraints.min_intersecting_angle_deg_other));
    cOther.setMinIntersectionMagnitude(cfg.constraints.min_intersection_magnitude);
    cOther.setMinIntersectionDistance(cfg.constraints.min_intersection_distance);

    EntityNetworkConstraintsMatrix<Constraints> constraintsMatrix;

    // Self constraints for each family
    for (std::size_t i = 0; i < cfg.families.size(); ++i)
    {
        const Id id(static_cast<int>(i + 1));
        constraintsMatrix.addConstraints(cSelf, IdPair(id, id));
    }

    // Cross constraints for each ordered pair i != j
    for (std::size_t i = 0; i < cfg.families.size(); ++i)
    {
        for (std::size_t j = 0; j < cfg.families.size(); ++j)
        {
            if (i == j)
                continue;
            constraintsMatrix.addConstraints(cOther, IdPair(Id(i + 1), Id(j + 1)));
        }
    }

    // initialize empty fracture families
    MultiGeometryEntitySet<Disk, Quad> entitySets;

    std::cout << "\n --- Start entity sampling ---\n"
              << std::endl;

    // Keep track of accepted entities
    std::vector<std::shared_ptr<Disk>> acceptedDisks;
    std::vector<std::shared_ptr<Quad>> acceptedQuads;

    while (!status.finished())
    {
        Id id;
        auto geom = multiSampler(id);

        if (status.finished(id))
        {
            status.increaseRejectedCounter("set finished");
            continue;
        }

        // constraints vs entities across all sets
        if (const auto res = constraintsMatrix.evaluate(entitySets, geom, id); res.violationDetected())
        {
            status.increaseRejectedCounter(res.violationLabel());
            continue;
        }

        entitySets.addEntity(geom, id);

        // store for CSV output (and metadata)
        if (auto d = std::dynamic_pointer_cast<Disk>(geom))
        {
            acceptedDisks.push_back(d);
        }
        else if (auto q = std::dynamic_pointer_cast<Quad>(geom))
        {
            acceptedQuads.push_back(q);
        }
        else
        {
            status.increaseRejectedCounter("unknown geometry type");
            continue;
        }

        status.increaseCounter(id);
        status.print();
    }

    std::cout << "\n --- Finished entity sampling ---\n"
              << std::endl;

    // Native Frackit CSV format:
    /*
    std::ofstream out(cfg.output.disks_csv);
    out << "family,cx,cy,cz,majx,majy,majz,minx,miny,minz,majR,minR,nx,ny,nz\n";

    for (std::size_t i = 0; i < entitySets.size(); ++i)
    {
        const int family = families[i].familyIndex;
        for (const auto &disk : entitySets[i])
        {
            const auto c = disk.center();
            const auto n = disk.normal();
            const auto maj = disk.majorAxis();
            const auto min = disk.minorAxis();
            out << family << ","
                << c.x() << "," << c.y() << "," << c.z() << ","
                << maj.x() << "," << maj.y() << "," << maj.z() << ","
                << min.x() << "," << min.y() << "," << min.z() << ","
                << disk.majorAxisLength() * 0.5 << ","
                << disk.minorAxisLength() * 0.5 << ","
                << n.x() << "," << n.y() << "," << n.z()
                << "\n";
        }
    }
    */

    // PorePy CSV format for disks:
    std::ofstream out_disks(cfg.output.disks_csv);

    // Start with the domain.
    // out << "DOMAIN_XMIN, DOMAIN_YMIN, DOMAIN_ZMIN, DOMAIN_XMAX, DOMAIN_YMAX, DOMAIN_ZMAX\n";
    out_disks << cfg.domain.xmin << "," << cfg.domain.ymin << "," << cfg.domain.zmin << ","
              << cfg.domain.xmax << "," << cfg.domain.ymax << "," << cfg.domain.zmax << "\n";

    // Now add the disks in the following format.
    // out << "CENTER_X, CENTER_X, CENTER_Y, CENTER_Z, MAJOR_AXIS, MINOR_AXIS, MAJOR_AXIS_ANGLE, STRIKE_ANGLE, DIP_ANGLE

    for (const auto &disk : acceptedDisks)
    {
        const auto c = disk->center();
        const auto n = disk->normal();
        const auto maj = disk->majorAxis();

        Vec3 n3{n.x(), n.y(), n.z()};
        Vec3 a3{maj.x(), maj.y(), maj.z()};

        auto [major_axis_angle, strike_angle, dip_angle] = toPorePyAngles(a3, n3);

        out_disks << c.x() << "," << c.y() << "," << c.z() << ","
                  << disk->majorAxisLength() * 0.5 << ","
                  << disk->minorAxisLength() * 0.5 << ","
                  << major_axis_angle << "," << strike_angle << "," << dip_angle << "\n";
    }

    // PorePy CSV format for quadrilaterals.
    std::ofstream out_quads(cfg.output.quads_csv);

    // Start with the domain.
    // out << "DOMAIN_XMIN, DOMAIN_YMIN, DOMAIN_ZMIN, DOMAIN_XMAX, DOMAIN_YMAX, DOMAIN_ZMAX\n";
    out_quads << cfg.domain.xmin << "," << cfg.domain.ymin << "," << cfg.domain.zmin << ","
              << cfg.domain.xmax << "," << cfg.domain.ymax << "," << cfg.domain.zmax << "\n";

    // Now add the disks in the following format (one line per fracture).
    // out << P0_X, P0_Y, P0_Z, ..., PN_X, PN_Y, PN_Z

    for (const auto &quad : acceptedQuads)
    {
        const auto p0 = quad->corner(0);
        const auto p1 = quad->corner(1);
        const auto p2 = quad->corner(2);
        const auto p3 = quad->corner(3);

        out_quads << p0.x() << "," << p0.y() << "," << p0.z() << ","
                  << p1.x() << "," << p1.y() << "," << p1.z() << ","
                  << p2.x() << "," << p2.y() << "," << p2.z() << ","
                  << p3.x() << "," << p3.y() << "," << p3.z() << "\n";
    }

    // Also write some metadata about the families to a separate file (optional, but can be useful for later reference)
    std::ofstream meta_out(cfg.output.families_csv);
    meta_out << "family_id,type,target_num,sampled_num\n";
    for (std::size_t i = 0; i < cfg.families.size(); ++i)
    {
        const auto &fam = cfg.families[i];
        meta_out << (i + 1) << "," << fam.type << "," << fam.target_num << "," << status.getCount(Id(static_cast<int>(i + 1))) << "\n";
    }

    // Write to monitor how many entities were sampled.
    for (std::size_t i = 0; i < cfg.families.size(); ++i)
        std::cout << "[monitor] Sampled entities (Family " << (i + 1)
                  << " with type " << cfg.families[i].type << "): "
                  << status.getCount(Id(static_cast<int>(i + 1))) << "\n";

    // We can now create an entity network from the two sets
    EntityNetworkBuilder<ctype> builder;
    entitySets.exportEntitySets(builder);
    const auto network = builder.build();

    // This can be written out in Gmsh (.geo) format to be
    // meshed by a two-dimensional surface mesh
    std::cout << "\n --- Writing .geo file ---\n"
              << std::endl;
    GmshWriter writer(network);
    // Set some mesh-size, but actually not important as the csv file is the main output here.
    // Export the .geo file for visualization purposes.
    writer.setMeshSize(GmshWriter::GeometryTag::entity, 0.1); // mesh size at entities
    writer.write("network");                                  // filename of the .geo files (will add extension .geo automatically)
    std::cout << "\n --- Finished writing .geo file ---\n"
              << std::endl;

    return 0;
}
