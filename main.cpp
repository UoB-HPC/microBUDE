
#include <cstddef>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <vector>
#include <array>
#include <chrono>
#include <string>

#ifndef PPWI
#define PPWI (64)
#endif

#ifndef ITERS
#define ITERS 8
#endif

#ifndef DECK_PATH
#define DECK_PATH  "bm1"
#endif


#define DIFF_TOLERANCE_PCT 0.025f
#define ENERGY_ENTRIES 8
#define FILE_LIGAND "/ligand.in"
#define FILE_PROTEIN "/protein.in"
#define FILE_FORCEFIELD "/forcefield.in"
#define FILE_POSES "/poses.in"
#define FILE_REF_ENERGIES "/ref_energies.out"

#define ZERO 0.0f
#define QUARTER 0.25f
#define HALF 0.5f
#define ONE 1.0f
#define TWO 2.0f
#define FOUR 4.0f
#define CNSTNT 45.0f

// Energy evaluation parameters
#define HBTYPE_F 70
#define HBTYPE_E 69
#define HARDNESS 38.0f
#define NPNPDIST 5.5f
#define NPPDIST 1.0f


struct __attribute__((__packed__)) Atom {
    float x, y, z;
    int32_t type;
};

struct __attribute__((__packed__)) FFParams {
    int32_t hbtype;
    float radius;
    float hphb;
    float elsc;
};

extern "C" void
fasten_main(size_t group, size_t ntypes, size_t nposes, size_t natlig, size_t natpro,        //
            const Atom *protein, const Atom *ligand,                                         //
            const float *transforms_0, const float *transforms_1, const float *transforms_2, //
            const float *transforms_3, const float *transforms_4, const float *transforms_5, //
            const FFParams *forcefield, float *energies                                      //
) {

    float transform[3][4][PPWI];
    float etot[PPWI];

#pragma omp simd
    for (int l = 0; l < PPWI; l++) {
        int ix = group * PPWI + l;

        // Compute transformation matrix
        const float sx = sinf(transforms_0[ix]);
        const float cx = cosf(transforms_0[ix]);
        const float sy = sinf(transforms_1[ix]);
        const float cy = cosf(transforms_1[ix]);
        const float sz = sinf(transforms_2[ix]);
        const float cz = cosf(transforms_2[ix]);

        transform[0][0][l] = cy * cz;
        transform[0][1][l] = sx * sy * cz - cx * sz;
        transform[0][2][l] = cx * sy * cz + sx * sz;
        transform[0][3][l] = transforms_3[ix];
        transform[1][0][l] = cy * sz;
        transform[1][1][l] = sx * sy * sz + cx * cz;
        transform[1][2][l] = cx * sy * sz - sx * cz;
        transform[1][3][l] = transforms_4[ix];
        transform[2][0][l] = -sy;
        transform[2][1][l] = sx * cy;
        transform[2][2][l] = cx * cy;
        transform[2][3][l] = transforms_5[ix];

        etot[l] = 0.f;
    }

    // Loop over ligand atoms
    for (int il = 0; il < natlig; il++) {
        // Load ligand atom data
        const Atom l_atom = ligand[il];
        const FFParams l_params = forcefield[l_atom.type];
        const int lhphb_ltz = l_params.hphb < 0.f;
        const int lhphb_gtz = l_params.hphb > 0.f;

        // Transform ligand atom
        float lpos_x[PPWI], lpos_y[PPWI], lpos_z[PPWI];

#pragma omp simd
        for (int l = 0; l < PPWI; l++) {
            lpos_x[l] = transform[0][3][l] + l_atom.x * transform[0][0][l] +
                        l_atom.y * transform[0][1][l] +
                        l_atom.z * transform[0][2][l];
            lpos_y[l] = transform[1][3][l] + l_atom.x * transform[1][0][l] +
                        l_atom.y * transform[1][1][l] +
                        l_atom.z * transform[1][2][l];
            lpos_z[l] = transform[2][3][l] + l_atom.x * transform[2][0][l] +
                        l_atom.y * transform[2][1][l] +
                        l_atom.z * transform[2][2][l];
        }
        // Loop over protein atoms
        for (int ip = 0; ip < natpro; ip++) {
            // Load protein atom data
            const Atom p_atom = protein[ip];
            const FFParams p_params = forcefield[p_atom.type];

            const float radij = p_params.radius + l_params.radius;
            const float r_radij = ONE / radij;

            const float elcdst = (p_params.hbtype == HBTYPE_F && l_params.hbtype == HBTYPE_F) ? FOUR
                                                                                              : TWO;
            const float elcdst1 = (p_params.hbtype == HBTYPE_F && l_params.hbtype == HBTYPE_F)
                                  ? QUARTER : HALF;
            const int type_E = ((p_params.hbtype == HBTYPE_E || l_params.hbtype == HBTYPE_E));

            const int phphb_ltz = p_params.hphb < 0.f;
            const int phphb_gtz = p_params.hphb > 0.f;
            const int phphb_nz = p_params.hphb != 0.f;
            const float p_hphb = p_params.hphb * (phphb_ltz && lhphb_gtz ? -ONE : ONE);
            const float l_hphb = l_params.hphb * (phphb_gtz && lhphb_ltz ? -ONE : ONE);
            const float distdslv = (phphb_ltz ? (lhphb_ltz ? NPNPDIST : NPPDIST) : (lhphb_ltz
                                                                                    ? NPPDIST
                                                                                    : -FLT_MAX));
            const float r_distdslv = ONE / distdslv;

            const float chrg_init = l_params.elsc * p_params.elsc;
            const float dslv_init = p_hphb + l_hphb;

#pragma omp simd
            for (int l = 0; l < PPWI; l++) {
                // Calculate distance between atoms
                const float x = lpos_x[l] - p_atom.x;
                const float y = lpos_y[l] - p_atom.y;
                const float z = lpos_z[l] - p_atom.z;
                const float distij = sqrtf(x * x + y * y + z * z);

                // Calculate the sum of the sphere radii
                const float distbb = distij - radij;

                const int zone1 = (distbb < ZERO);

                // Calculate steric energy
                etot[l] += (ONE - (distij * r_radij)) * (zone1 ? TWO * HARDNESS : 0.f);

                // Calculate formal and dipole charge interactions
                float chrg_e = chrg_init * ((zone1 ? ONE : (ONE - distbb * elcdst1)) *
                                            (distbb < elcdst ? ONE : ZERO));
                float neg_chrg_e = -fabsf(chrg_e);
                chrg_e = type_E ? neg_chrg_e : chrg_e;
                etot[l] += chrg_e * CNSTNT;

                // Calculate the two cases for Nonpolar-Polar repulsive interactions
                float coeff = (ONE - (distbb * r_distdslv));
                float dslv_e = dslv_init * ((distbb < distdslv && phphb_nz) ? ONE : 0.f);
                dslv_e *= (zone1 ? ONE : coeff);
                etot[l] += dslv_e;
            }
        }
    }

    // Write result
#pragma omp simd
    for (int l = 0; l < PPWI; l++) {
        energies[group * PPWI + l] = etot[l] * HALF;
    }
}

template<typename T>
std::vector<T> readNStruct(const std::string &path) {
    FILE *fp = fopen(path.data(), "rb");
    if (!fp) {
        fprintf(stderr, "Bad file: %s\n", path.data());
        std::exit(1);
    }
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);
    rewind(fp);

    std::vector<T> xs(len / sizeof(T));
    if (fread(xs.data(), sizeof(T), len / sizeof(T), fp) != len / sizeof(T)) {
        fprintf(stderr, "Failed to read file\n");
        std::exit(1);
    }
    fclose(fp);
    return xs;
}

int main() {

    const std::string deckDir = DECK_PATH;

    auto ligand = readNStruct<Atom>(deckDir + FILE_LIGAND);
    auto protein = readNStruct<Atom>(deckDir + FILE_PROTEIN);
    auto forcefield = readNStruct<FFParams>(deckDir + FILE_FORCEFIELD);

    auto poseData = readNStruct<float>(deckDir + FILE_POSES);
    if (poseData.size() % 6 != 0) {
        fprintf(stderr, "Pose size (%zu) not divisible my 6!", poseData.size());
        std::exit(1);
    }
    auto maxPoses = poseData.size() / 6;
    auto nposes = maxPoses;

    if (nposes > maxPoses) {
        fprintf(stderr, "Requested poses (%zu) exceeded max poses (%zu) for deck\n", nposes,
                maxPoses);
        std::exit(1);
    }
    if (nposes < PPWI) {
        fprintf(stderr, "pose count %zu <= %d (PPWI)\n", nposes, PPWI);
        std::exit(1);
    }

    if (nposes % PPWI != 0) {
        fprintf(stderr, "pose count %zu %% %d (PPWI) != 0", nposes, PPWI);
        std::exit(1);
    }

    std::array<std::vector<float>, 6> poses;
    for (size_t i = 0; i < 6; ++i) {
        poses[i].resize(nposes);
        std::copy(std::next(poseData.cbegin(), int(i * maxPoses)),
                  std::next(poseData.cbegin(), int(i * maxPoses + nposes)),
                  poses[i].begin());
    }

    std::array<float *, 6> data_poses{};
    auto data_protein = static_cast<Atom *>(std::malloc(sizeof(Atom) * protein.size()));
    auto data_ligand = static_cast<Atom *>(std::malloc(sizeof(Atom) * ligand.size()));
    auto data_forcefield = static_cast<FFParams *>(std::malloc(
            sizeof(FFParams) * forcefield.size()));
    auto energies = static_cast<float *>(std::malloc(sizeof(float) * nposes));

    for (auto i = 0; i < 6; i++)
        data_poses[i] = static_cast<float *>(std::malloc(sizeof(float) * nposes));

#pragma omp parallel
    {
        for (auto i = 0; i < 6; i++) {
#pragma omp for nowait
            for (auto j = 0; j < nposes; j++)
                data_poses[i][j] = poses[i][j];
        }
#pragma omp for nowait
        for (auto i = 0; i < nposes; i++)
            energies[i] = 0.f;

#pragma omp for nowait
        for (auto i = 0; i < protein.size(); i++)
            data_protein[i] = protein[i];

#pragma omp for nowait
        for (auto i = 0; i < ligand.size(); i++)
            data_ligand[i] = ligand[i];

#pragma omp for nowait
        for (auto i = 0; i < forcefield.size(); i++)
            data_forcefield[i] = forcefield[i];
    }

    auto poses_0 = data_poses[0];
    auto poses_1 = data_poses[1];
    auto poses_2 = data_poses[2];
    auto poses_3 = data_poses[3];
    auto poses_4 = data_poses[4];
    auto poses_5 = data_poses[5];

    float minRuntimeS = std::numeric_limits<float>::max();

    printf("PPWI=%d\n", PPWI);
    printf("Deck=%s\n", deckDir.data());

    for (int iter = 0; iter < ITERS; ++iter) {
        auto start = std::chrono::high_resolution_clock::now();


#pragma omp parallel for
        for (size_t group = 0; group < (nposes / PPWI); group++) {
            fasten_main(group, forcefield.size(), nposes, ligand.size(),
                        protein.size(),                //
                        data_protein, data_ligand,                                      //
                        poses_0, poses_1, poses_2, poses_3, poses_4, poses_5, //
                        data_forcefield, energies);
        }

        auto end = std::chrono::high_resolution_clock::now();


        std::chrono::duration<float, std::milli> elapsedMs = end - start;

        printf("[%d] = %f ms\n", iter, elapsedMs.count());
        minRuntimeS = fminf(minRuntimeS, elapsedMs.count());
    }

    printf("Min = %f ms\n", minRuntimeS);


    auto refEnergiesData = readNStruct<char>(deckDir + FILE_REF_ENERGIES);
    std::string refEnergies(refEnergiesData.begin(), refEnergiesData.end());
    size_t pos = 0, prev = 0, entry = 0;
    while ((pos = refEnergies.find("\n", prev)) != std::string::npos) {
        auto ref = std::stof(refEnergies.substr(prev, pos - prev));
        auto diffPct = (fabsf(ref - energies[entry]) / ref) * 100.f;
        if (entry < ENERGY_ENTRIES)
            printf("[%zu] %f diff=(%f)\n", entry, energies[entry], diffPct);
        if (diffPct > DIFF_TOLERANCE_PCT && ref >= 1.f && energies[entry] >= 1.f)
            fprintf(stderr, "[%zu] Invalid: %f diff=(%f)\n", entry, energies[entry], diffPct);
        prev = pos + 1;
        entry++;
    }


    return EXIT_SUCCESS;
}
