/*
 * =============================================================================
 * mirna.h
 *
 * Struct used for storing miRNA information.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_MIRNA_H
#define ISOMIR_SEA_MIRNA_H

// ============================================================================
// Prerequisites
// ============================================================================

unsigned const END5P = 0;
unsigned const END3P = 1;
unsigned const ENDUNK = 2;

unsigned MAX_MIR_SIZE = 24; //TODO verify this value

bool const MATURE = false;
bool const STAR = true;

typedef struct t_mir_info_cell
{
    std::string info{}; // full information
    std::string mir_id{}; // mirgenedb id
    std::string mimat{}; // mirbase id
    std::string mimat_prx{}; // mirbase id of precursor
    unsigned end5p3p{ENDUNK};
    bool mature_star{MATURE};
    int index{-1};
    // GFF3 fields
    std::string ref{}; // chromosome
    int start_gen{-1};
    int end_gen{-1};
    char strand{'?'}; // false: - | true: +
    int count_overwrite{0};

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive( info, mir_id, mimat, mimat_prx, end5p3p, mature_star, index, ref, start_gen, end_gen, strand, count_overwrite );
    }
} t_mir_info_cell;

#endif //ISOMIR_SEA_MIRNA_H
