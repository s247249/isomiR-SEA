/*
 * =============================================================================
 * tag_mirna.h
 *
 * Structs used for storing the tag-miRNA matching fields.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_TAG_MIRNA_H
#define ISOMIR_SEA_TAG_MIRNA_H

unsigned const START_CENTRAL_SITE = 3;
unsigned const END_CENTRAL_SITE = 17;
unsigned const SIZE_CENTRAL_SITE = 11;
unsigned const START_SUPP_SITE = 12;
unsigned const END_SUPP_SITE = 17;
unsigned const SIZE_SUPP_SITE = 3;
unsigned const START_COMP_SITE = 11;
unsigned const END_COMP_SITE = 21;
unsigned const SIZE_COMP_SITE = 8;

unsigned const INDEL_PENALTY = 1;
unsigned const MISMATCH_UNALIGN_PENALTY = 2;
unsigned const MISMATCH_ALIGN_PENALTY = 1;
unsigned const MISMATCH_ALIGN_PENALTY_SEED = 3;

std::vector<std::pair<std::string, std::string> > tmf = // tagMirFields
        {
                // TAG
                {"tag_id", "TI"},                  //  0
                {"tag_seq", "TS"},
                {"tag_qual", "TQ"},
                {"tag_count", "TC"},
                {"tag_start_gen", "TSG"},
                {"tag_end_gen", "TEG"},             //  5
                // ORGANISM
                {"org_code", "ORG"},
                // MIR SEQ
                {"mir_id", "MI"},
                {"mir_seq", "MS"},
                // ALIGN
                {"align_mark", "AM"},
                {"align_iupac_tag", "IT"},           // 10
                {"align_iupac_mir", "IM"},
                {"align_cigar", "CI"},
                {"align_length", "AL"},
                {"mir_tag_size_diff", "SD"},
                {"iso_exact", "IEX"},               // 15
                {"iso5p", "I5P"},
                {"iso_m_snp", "IMS"},
                {"iso_snp", "ISN"},
                {"iso3p", "I3P"},
                {"mismatch_in_seed", "INS"},         // 20
                {"off_site", "IOS"},
                {"suppl_site", "ISS"},
                {"compens_site", "IPS"},
                {"central_site", "ICS"},
                {"a2i_in_seed", "AIS"},              // 25
                {"a2i_out_seed", "AIO"},
                // Align score difference and sizes mir
                {"num_of_mir_same_seq", "MI4S"},
                {"score_diff_mir0_mir1", "MSD"},
                // MIR INFO
                {"mir_info_index", "MII"},
                {"mir_info_info", "MIN"},            // 30
                {"mir_info_ref", "MRF"},
                {"mir_info_start_gen", "MSG"},
                {"mir_info_end_gen", "MEG"},
                {"mir_info_mimat", "MIM"},
                {"mir_info_strand", "MSR"},          // 35
                // PRXMIR
                {"prx_index", "PII"},
                {"prx_info", "PIN"},
                {"prx_ref", "PRF"},
                {"prx_seq", "PS"},
                {"prx_start_gen", "PSG"},            // 40
                {"prx_end_gen", "PEG"},
                {"prx_mimat", "PMI"},
                {"prx_strand", "PSR"},
                // Align score difference and sizes premir
                {"num_of_prx_same_seq", "PR4S"},
                {"prx_align_mark", "PAM"},      // 45
                {"iso5p_canonic", "IC5"},
                {"iso3p_canonic", "IC3"},
                {"score_diff_prx0_prx1", "PSD"},
                // INDEX
                {"mir_tag_out_id", "MTOID"},          // 49
                {"mir_matched_id", "MMI"}, // 50 TODO REMOVE after checking if equal to MI This should be the same as MI (7)
                {"mir_info_id", "MIID"}, //   51 TODO REMOVE This should be the same as MII (29)
                {"premir_id", "PID"}, //      52 TODO REMOVE This should be the same as PII (36)
                {"align_iupac_prx", "IP"},
                {"single_multi_discarded", "SMD"}
        };

typedef struct t_tag_mir_matched_v
{
    t_mir_matched_v mir_matched;
    t_tag_cell * tag;
}t_tag_mir_matched_v;

typedef std::vector<t_tag_mir_matched_v> t_tag_mir_matched_vv; // used to save the tag-mir associations

typedef std::variant<int, unsigned, char, std::string, t_seq> t_variant;
typedef std::map<std::string, std::pair<std::string, t_variant> > t_map_str_any;
//typedef std::vector<t_map_str_any> t_vect_map_str_any;

#endif //ISOMIR_SEA_TAG_MIRNA_H
