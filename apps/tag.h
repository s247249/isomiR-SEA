/*
 * =============================================================================
 * tag.h
 *
 * Struct used for storing tag information.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_TAG_H
#define ISOMIR_SEA_TAG_H

typedef struct t_tag_cell
{
    t_seq seq{};
    std::string qual{};
    unsigned count{};
    unsigned index{};
    int mir_matched_id{-1};
    int prxmir_matched_id{-1};
    std::string reason_discard{};
}t_tag_cell;

typedef std::vector<t_tag_cell > t_tag_v;

typedef struct t_tag
{
    unsigned tag_count{};
    unsigned read_count{}; // counter for check how many reads have been analized
    unsigned invalid_tag_count{};
    unsigned invalid_read_count{};
    t_tag_v tag_v;
    t_tag_v discarded_tag_v;
}t_tag;

#endif //ISOMIR_SEA_TAG_H
