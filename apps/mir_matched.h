/*
 * =============================================================================
 * mir_matched.h
 *
 * Struct used for storing the read mapping results.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_MIR_MATCHED_H
#define ISOMIR_SEA_MIR_MATCHED_H

typedef struct t_mir_matched_cell // type used for store the read mapping results
{
    std::map<std::string, t_prxmir_matched_v> mir_prxmir;
    const t_mir_cell * mir{NULL};
    int mismatch_in_seed{};
    int align_score{};
    int align_begin_mir;
    int align_begin_tag;
    int align_size;
    int mir_tag_size_diff;
    // Boolean variable used to create more accurate miRNA expression level profile
    bool mir_exact{}; // canonical mature sequence (this can be extracted from the difference with the other cases)
    int iso5p{}; // 5p isoform, characterized by insertions or deletions in the 5’ end of miRNA
    int iso3p{}; // 3p isoform, characterized by insertions or deletions in the 3’ end of miRNA
    bool iso_snp{}; // single nucleotide polymorphism isoforms, presenting a mismatch with respect to miRNA sequence
    bool iso_msnp{}; // multiple nucleotide polymorphism isoform, presenting more than a mismatch with respect to miRNA sequence
    bool off_site{};
    bool suppl_site{};
    bool compens_site{};
    bool central_site{}; // for now is the inclusion of central site 3-13, 4-14 and 5-15
    std::string cigar;
    std::string align_iupac_mir;
    std::string align_iupac_tag;
    bool a2i_seed{};
    bool a2i_out_seed{};
} t_mir_matched_cell;

typedef std::vector<t_mir_matched_cell> t_mir_matched_v;

#endif //ISOMIR_SEA_MIR_MATCHED_H
