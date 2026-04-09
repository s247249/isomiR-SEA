/*
 * =============================================================================
 * seed_mirna.h
 *
 * Struct used for storing seed and miRNA sequences.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_SEED_MIRNA_H
#define ISOMIR_SEA_SEED_MIRNA_H

typedef struct t_mir_cell
{
    t_seq seq;
    unsigned index; // this index is automatically generated when the reference is acquired
    std::map<std::string, t_prx_mirna_v> org_prx_m; // full information grouped 4 organism

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive( seq, index, org_prx_m );
    }
} t_mir_cell;

typedef std::vector<t_mir_cell> t_mirna_v;

typedef struct t_seed_cell
{
    t_seq seq;
    unsigned index; // automatically generated when the reference is acquired
    t_mirna_v mir_seqs;

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive( seq, index, mir_seqs );
    }
} t_seed_cell;

typedef std::vector<t_seed_cell > t_seed_v;

typedef struct t_seed
{
    int bph; // begin position horizontal/miRNA/database
    int bpv; // begin position vertical/tag/query
    int eph; // end position horizontal/miRNA/database
    int epv; // end position vertical/tag/query
} t_seed;

#endif //ISOMIR_SEA_SEED_MIRNA_H
