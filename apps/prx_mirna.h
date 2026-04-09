/*
 * =============================================================================
 * prx_mirna.h
 *
 * Struct used for storing prxmiRNA information.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_PRX_MIRNA_H
#define ISOMIR_SEA_PRX_MIRNA_H

using namespace seqan3;

int const ORGCODELENGTH = 3;
std::string const ALL_ORGANISM = "all";

typedef std::vector<rna15> t_seq;
typedef std::vector<std::vector<phred42>> t_qual;
typedef std::vector<std::string> t_vect_str;
typedef std::map<std::string, bool> t_map_str_bool;
typedef std::pair<std::string, bool> t_pair_str_bool;

int const MIRNA_SIZE = 22;
int const FLANK_SIZE = 30;
int const LOOP_SIZE = 30;

t_seq mirna_N(MIRNA_SIZE, 'N'_rna15);
t_seq flank_N(FLANK_SIZE, 'N'_rna15);
t_seq loop_N(LOOP_SIZE, 'N'_rna15);

typedef struct t_oth_loc // other location of this prxmir in the genomes
{
    std::string name{}; // mir_id or mimat
    int start_gen_pre{-1};
    int end_gen_pre{-1};
    std::string ref{}; // chromosome

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive( name, start_gen_pre, end_gen_pre, ref );
    }
} t_oth_loc;

typedef struct t_prx_mirna_cell
{
    t_seq prx_seq{}; // can be a primir or a premir sequence
    t_seq tmp5p{};
    t_seq tmp3p{};
    std::string info{}; // full information
    std::string mir_id{}; // mirgenedb id
    std::string mimat{}; // mirbase id
    int index{-1};
    int index5p{-1};
    int index3p{-1};
    // GFF3 fields
    std::string ref{}; // chromosome
    int start_gen_pri{-1};
    int end_gen_pri{-1};
    int start_gen_pre{-1};
    int end_gen_pre{-1};
    char strand{'?'}; // false = - | true = +
    std::vector<t_oth_loc> other_locations;
    // 5p and 3p
    t_mir_info_cell mirna5p;
    t_mir_info_cell mirna3p;
    t_mir_info_cell mirna_tmp;

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive( prx_seq, tmp5p, tmp3p, info, mir_id, mimat, index, index5p, index3p, ref, start_gen_pri, end_gen_pri, start_gen_pre, end_gen_pre, strand, other_locations, mirna5p, mirna3p, mirna_tmp );
    }
} t_prx_mirna_cell;

typedef std::vector<t_prx_mirna_cell > t_prx_mirna_v;
typedef std::map<std::string, t_prx_mirna_cell> t_prx_mir_m; // map<prxmir_id, t_prx_mirna_cell>
typedef std::map<std::string, t_prx_mir_m> t_org_prxmir_m; // map<organism, t_prx_mir_m>

typedef std::pair<std::string, t_prx_mirna_cell> t_prx_mir_p; // pair<organism, t_prx_mirna_cell> used to debug-print prxmir cells
typedef std::pair<std::string, t_prx_mir_m> t_org_prxmir_p; // pair<organism, t_prx_mir_m> used to debug-print prxmir cells
typedef std::pair<std::string, t_prx_mirna_v> t_org_prx_p; // pair<organism, t_prx_mirna_v> used to debug-print prxmir cells starting from seeds

#endif //ISOMIR_SEA_PRX_MIRNA_H
