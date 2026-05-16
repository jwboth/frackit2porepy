#include <random>
#include <iostream>
#include <fstream>
#include <string>

#include <frackit/common/id.hh>
#include <frackit/io/gmshwriter.hh>

#include <frackit/sampling/makeuniformpointsampler.hh>
#include <frackit/sampling/disksampler.hh>
#include <frackit/sampling/quadrilateralsampler.hh>
#include <frackit/sampling/status.hh>

#include <frackit/geometry/quadrilateral.hh>
#include <frackit/geometry/box.hh>

#include <frackit/entitynetwork/constraints.hh>
#include <frackit/entitynetwork/networkbuilder.hh>

#include "toml.hpp" // Third-party TOML library

namespace
{

    struct NormalParams
    {
        double mean = 0.0;
        double stddev = 0.0;
    };

    struct DomainConfig
    {
        double xmin = 0, ymin = 0, zmin = 0, xmax = 100, ymax = 100, zmax = 100;
    };

    struct SamplerConfig
    {
        NormalParams major_axis_length{50.0, 10.0};
        NormalParams minor_axis_length{50.0, 5.0};
        NormalParams rot_x_deg{0.0, 7.5};
        NormalParams rot_y_deg{0.0, 7.5};
        NormalParams rot_z_deg{0.0, 7.5};
    };

    struct ConstraintsConfig
    {
        double min_distance = 0.05;
        double min_intersecting_angle_deg_self = 30.0;
        double min_intersecting_angle_deg_other = 40.0;
        double min_intersection_magnitude = 0.05;
        double min_intersection_distance = 0.05;
    };

    struct SamplingConfig
    {
        int target_num_disks = 5;
    };

    struct OutputConfig
    {
        std::string disks_csv = "disks.csv";
    };

    struct Config
    {
        DomainConfig domain{};
        SamplerConfig sampler1{};
        SamplerConfig sampler2{SamplerConfig{
            NormalParams{50.0, 10.0}, NormalParams{50.0, 5.0},
            NormalParams{45.0, 7.5}, NormalParams{60.0, 7.5}, NormalParams{90.0, 7.5}}};
        ConstraintsConfig constraints{};
        SamplingConfig sampling{};
        OutputConfig output{};
    };

    static void printConfig(const Config &cfg)
    {
        std::cout << "\n[config] domain: [" << cfg.domain.xmin << ", " << cfg.domain.ymin << ", " << cfg.domain.zmin
                  << "] -> [" << cfg.domain.xmax << ", " << cfg.domain.ymax << ", " << cfg.domain.zmax << "]\n";
        auto ps = [](const char *name, const SamplerConfig &s)
        {
            std::cout << "[config] " << name << ".major_axis_length: mean=" << s.major_axis_length.mean << " stddev=" << s.major_axis_length.stddev << "\n";
            std::cout << "[config] " << name << ".minor_axis_length: mean=" << s.minor_axis_length.mean << " stddev=" << s.minor_axis_length.stddev << "\n";
            std::cout << "[config] " << name << ".rot_x_deg: mean=" << s.rot_x_deg.mean << " stddev=" << s.rot_x_deg.stddev << "\n";
            std::cout << "[config] " << name << ".rot_y_deg: mean=" << s.rot_y_deg.mean << " stddev=" << s.rot_y_deg.stddev << "\n";
            std::cout << "[config] " << name << ".rot_z_deg: mean=" << s.rot_z_deg.mean << " stddev=" << s.rot_z_deg.stddev << "\n";
        };
        ps("sampler1", cfg.sampler1);
        ps("sampler2", cfg.sampler2);

        std::cout << "[config] constraints.min_distance=" << cfg.constraints.min_distance << "\n";
        std::cout << "[config] constraints.min_intersecting_angle_deg_self=" << cfg.constraints.min_intersecting_angle_deg_self << "\n";
        std::cout << "[config] constraints.min_intersecting_angle_deg_other=" << cfg.constraints.min_intersecting_angle_deg_other << "\n";
        std::cout << "[config] constraints.min_intersection_magnitude=" << cfg.constraints.min_intersection_magnitude << "\n";
        std::cout << "[config] constraints.min_intersection_distance=" << cfg.constraints.min_intersection_distance << "\n";
        std::cout << "[config] sampling.target_num_disks=" << cfg.sampling.target_num_disks << "\n";
        std::cout << "[config] output.disks_csv=\"" << cfg.output.disks_csv << "\"\n";
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

            auto loadSampler = [](const toml::table &root, const char *key, SamplerConfig &s)
            {
                const auto *st = root[key].as_table();
                if (!st)
                    return;

                if (const auto *maj = (*st)["major_axis_length"].as_table())
                {
                    s.major_axis_length.mean = (*maj)["mean"].value_or(s.major_axis_length.mean);
                    s.major_axis_length.stddev = (*maj)["stddev"].value_or(s.major_axis_length.stddev);
                }
                if (const auto *min = (*st)["minor_axis_length"].as_table())
                {
                    s.minor_axis_length.mean = (*min)["mean"].value_or(s.minor_axis_length.mean);
                    s.minor_axis_length.stddev = (*min)["stddev"].value_or(s.minor_axis_length.stddev);
                }

                if (const auto *rot = (*st)["rotation_deg"].as_table())
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

            loadSampler(tbl, "sampler1", cfg.sampler1);
            loadSampler(tbl, "sampler2", cfg.sampler2);

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

            if (const auto *s = tbl["sampling"].as_table())
                cfg.sampling.target_num_disks = (*s)["target_num_disks"].value_or(cfg.sampling.target_num_disks);

            if (const auto *o = tbl["output"].as_table())
                cfg.output.disks_csv = (*o)["disks_csv"].value_or(cfg.output.disks_csv);
        }
        catch (const toml::parse_error &err)
        {
            std::cerr << "[config] TOML parse error in " << path << ": " << err.description() << "\n";
        }
    }

} // namespace

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

    // Internal geometry type for 3d quadrilaterals
    // using Quad = Quadrilateral<ctype, worldDimension>;
    using Disk = Disk<ctype>;

    // Define a domain (here: unit cube) in which the quadrilaterals should be created.
    // Boxes are created by providing xmin, ymin, zmin and xmax, ymax and zmax in constructor.
    Box<ctype> domain(cfg.domain.xmin, cfg.domain.ymin, cfg.domain.zmin,
                      cfg.domain.xmax, cfg.domain.ymax, cfg.domain.zmax);

    // We now create a sampler instance that uniformly samples points within this box.
    // These points will be used as the quadrilateral centers.
    auto pointSampler = makeUniformPointSampler(domain);

    // Sampler class for quadrilaterals. Per default, this uses uniform distributions
    // for all parameters defining the quadrilaterals. Quadrilateral samplers require
    // distributions for strike angle, dip angle, edge length (see class description
    // for more details).
    using NormalDistro = std::normal_distribution<ctype>;
    using UniformDistro = std::uniform_real_distribution<ctype>;

    // sampler for disks (orientation 1)
    DiskSampler diskSampler1(
        pointSampler,
        NormalDistro(cfg.sampler1.major_axis_length.mean, cfg.sampler1.major_axis_length.stddev),
        NormalDistro(cfg.sampler1.minor_axis_length.mean, cfg.sampler1.minor_axis_length.stddev),
        NormalDistro(toRadians(cfg.sampler1.rot_x_deg.mean), toRadians(cfg.sampler1.rot_x_deg.stddev)),
        NormalDistro(toRadians(cfg.sampler1.rot_y_deg.mean), toRadians(cfg.sampler1.rot_y_deg.stddev)),
        NormalDistro(toRadians(cfg.sampler1.rot_z_deg.mean), toRadians(cfg.sampler1.rot_z_deg.stddev)));

    // sampler for disks (orientation 1)
    DiskSampler diskSampler2(
        pointSampler,
        NormalDistro(cfg.sampler2.major_axis_length.mean, cfg.sampler2.major_axis_length.stddev),
        NormalDistro(cfg.sampler2.minor_axis_length.mean, cfg.sampler2.minor_axis_length.stddev),
        NormalDistro(toRadians(cfg.sampler2.rot_x_deg.mean), toRadians(cfg.sampler2.rot_x_deg.stddev)),
        NormalDistro(toRadians(cfg.sampler2.rot_y_deg.mean), toRadians(cfg.sampler2.rot_y_deg.stddev)),
        NormalDistro(toRadians(cfg.sampler2.rot_z_deg.mean), toRadians(cfg.sampler2.rot_z_deg.stddev)));

    // We want to enforce some constraints on the set of quadrilaterals.
    // In particular, for entities of the same set we want a minimum spacing
    // distance of 5cm, and the quadrilaterals must not intersect in angles
    // less than 30°. Moreover, if they intersect, we don't want intersection
    // edges whose length is smaller than 5cm, and, the intersection should not
    // be too close to the boundary of one of two intersecting quadrilaterals. Here: 5cm.
    EntityNetworkConstraints<ctype> constraintsOnSelf;
    constraintsOnSelf.setMinDistance(cfg.constraints.min_distance);
    constraintsOnSelf.setMinIntersectingAngle(toRadians(cfg.constraints.min_intersecting_angle_deg_self));
    constraintsOnSelf.setMinIntersectionMagnitude(cfg.constraints.min_intersection_magnitude);
    constraintsOnSelf.setMinIntersectionDistance(cfg.constraints.min_intersection_distance);

    // with respect to entities of the other set, we want to have larger intersection angles
    auto constraintsOnOther = constraintsOnSelf;
    constraintsOnOther.setMinIntersectingAngle(toRadians(cfg.constraints.min_intersecting_angle_deg_other));

    // container to store created entities
    std::vector<Disk> entitySet1;
    std::vector<Disk> entitySet2;

    // we give unique identifiers to both entity sets
    const Id idSet1(1);
    const Id idSet2(2);

    // use the status class to define when to stop sampling
    const auto target_num_disks = cfg.sampling.target_num_disks;
    SamplingStatus status;
    status.setTargetCount(idSet1, argc > 1 ? std::stoi(argv[1]) : target_num_disks); // per default, we want 5 entities in set 1
    status.setTargetCount(idSet2, argc > 1 ? std::stoi(argv[1]) : target_num_disks); // per default, we want 5 entities in set 2

    // start sampling into set 1 and keep alternating
    bool sampleIntoSet1 = true;

    std::cout << "\n --- Start entity sampling ---\n"
              << std::endl;
    while (!status.finished())
    {
        // sample a disk, alternating between sampler 1 and sampler 2
        auto disk = sampleIntoSet1 ? diskSampler1() : diskSampler2();
        auto &entitySet = sampleIntoSet1 ? entitySet1 : entitySet2;
        const auto &otherEntitySet = sampleIntoSet1 ? entitySet2 : entitySet1;

        // sample again if constraints w.r.t. other
        // entities of this set are not fulfilled
        if (!constraintsOnSelf.evaluate(entitySet, disk))
        {
            status.increaseRejectedCounter();
            continue;
        }

        // sample again if constraints w.r.t. other
        // entities of the other set are not fulfilled
        if (!constraintsOnOther.evaluate(otherEntitySet, disk))
        {
            status.increaseRejectedCounter();
            continue;
        }

        // the disk is admissible
        entitySet.push_back(disk);

        // print disk properties to terminal
        auto disk_family = sampleIntoSet1 ? 1 : 2;
        auto disk_center = disk.center();
        auto disk_major_axis = disk.majorAxis();
        auto disk_minor_axis = disk.minorAxis();
        auto disk_major_radius = disk.majorAxisLength() / 2.0;
        auto disk_minor_radius = disk.minorAxisLength() / 2.0;

        // tell the status we have a new entity in this set
        const auto &setId = sampleIntoSet1 ? idSet1 : idSet2;
        status.increaseCounter(setId);
        status.print();

        // sample into the other set the next time
        sampleIntoSet1 = !sampleIntoSet1;
    }
    std::cout << "\n --- Finished entity sampling ---\n"
              << std::endl;

    std::ofstream out(cfg.output.disks_csv);

    out << "family,cx,cy,cz,majx,majy,majz,minx,miny,minz,majR,minR,nx,ny,nz\n";

    auto dumpSet = [&](const std::vector<Disk> &set, int family)
    {
        for (const auto &disk : set)
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
    };

    dumpSet(entitySet1, 1);
    dumpSet(entitySet2, 2);

    // Write to monitor how many entities were sampled and how many were rejected due to constraint violations
    std::cout << "[monitor] Sampled entities (Set 1): " << entitySet1.size() << "\n";
    std::cout << "[monitor] Sampled entities (Set 2): " << entitySet2.size() << "\n";

    // We can now create an entity network from the two sets
    EntityNetworkBuilder<ctype> builder;
    builder.addEntities(entitySet1);
    builder.addEntities(entitySet2);

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
