/*
 * =============================================================================
 * mir_matched.h
 *
 * Struct used for storing the read mapping results (for prxmirs).
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_PRXMIR_MATCHED_H
#define ISOMIR_SEA_PRXMIR_MATCHED_H

typedef struct t_prxmir_matched_cell // type used for store the read mapping results
{
    const t_mir_info_cell * mir_info{NULL};
    const t_prx_mirna_cell * prxmir{NULL};
    std::string org{}; // TODO this can be removed
    int mir_info_id{-1};
    int prxmir_id{-1};
    int align_score{};
    bool iso3p_canonic{};
    bool iso5p_canonic{};
    int align_begin_prx{}; // used for tag alignment
    int align_begin_tag{}; // used for tag alignment
    int align_size{};
    std::string cigar;
    std::string align_iupac_prx;
    std::string align_iupac_tag;
}t_prxmir_matched_cell;

typedef std::vector<t_prxmir_matched_cell > t_prxmir_matched_v;
typedef std::pair<std::string, t_prxmir_matched_v> t_org_prxmir_matched_p;

#endif //ISOMIR_SEA_PRXMIR_MATCHED_H
