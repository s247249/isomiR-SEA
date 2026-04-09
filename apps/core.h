/*
 * =============================================================================
 * core.h
 *
 * Core of the tool (alignment algorithm).
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_CORE_H
#define ISOMIR_SEA_CORE_H

typedef seqan::Prefix<seqan::IupacString>::Type t_prefix;

enum ExtensionDirection
{
    EXTEND_NONE  = 0,
    EXTEND_LEFT  = 1,
    EXTEND_RIGHT = 2,
    EXTEND_BOTH  = 3
};

// ----------------------------------------------------------------------------
// Function from_rna15_to_iupac()
// ----------------------------------------------------------------------------

void from_rna15_to_iupac(seqan::IupacString & seq_s2, t_seq & seq_s3)
{
    for(rna15 c : seq_s3) {
        seqan::appendValue(seq_s2, c.to_char());
    }
}

template <typename t_cell>
bool compare_by_score(const t_cell & a, const t_cell & b)
{
    return a.align_score > b.align_score;
}

// ----------------------------------------------------------------------------
// Function extend_seed() adaptation of (old seqan) extendSeed <include/seqan/seeds/seeds_extension.h>
// ----------------------------------------------------------------------------

void extend_seed(t_seed & seed, t_seq const & database, t_seq const & query, ExtensionDirection direction)
{
    // Extension to the left
    if (direction == EXTEND_LEFT || direction == EXTEND_BOTH)
    {
        int posh = seed.bph;
        int posv = seed.bpv;
        while (posh >= 1 && posv >= 1 && database[posh - 1] == query[posv - 1])
        {
            --posh;
            --posv;
        }
        seed.bph = posh;
        seed.bpv = posv;
    }
    // Extension to the right
    if (direction == EXTEND_RIGHT || direction == EXTEND_BOTH)
    {
        int lenh = database.size();
        int lenv = query.size();
        int posh = seed.eph;
        int posv = seed.epv;
        while (posh < lenh && posv < lenv && database[posh] == query[posv])
        {
            ++posh;
            ++posv;
        }
        seed.eph = posh;
        seed.epv = posv;
    }
}

// ----------------------------------------------------------------------------
// Function seed_ext_both()
// ----------------------------------------------------------------------------

void seed_ext_both(int & bph, int & bpv, int & aln_size, t_seq const & mir_seq, t_seq const & tag_seq)
{
    if(bph != 0 && bpv != 0)
    {
        --bph;
        --bpv;
        ++aln_size;
    }
    if(mir_seq.size() > bph + aln_size && tag_seq.size() > bpv + aln_size)
    {
        ++aln_size;
    }
    t_seed seed{bph, bpv, bph+aln_size, bpv+aln_size};
    extend_seed(seed, mir_seq, tag_seq, EXTEND_BOTH);
    bph = seed.bph;
    bpv = seed.bpv;
    aln_size = seed.eph - bph;
}

// ----------------------------------------------------------------------------
// Function seed_ext_right()
// ----------------------------------------------------------------------------

void seed_ext_rigth(int & bph, int & bpv, int & aln_size, t_seq const & mir_seq, t_seq const & tag_seq)
{
    if(mir_seq.size() > bph + aln_size && tag_seq.size() > bpv + aln_size)
    {
        ++aln_size;
    }
    t_seed seed{bph, bpv, bph+aln_size, bpv+aln_size};
    extend_seed(seed, mir_seq, tag_seq, EXTEND_RIGHT);
    aln_size = seed.eph - bph;
}

// ----------------------------------------------------------------------------
// Function select_prxmir_ext()
// ----------------------------------------------------------------------------

void select_prxmir_ext(t_mir_matched_cell & mir_match, t_prxmir_matched_cell & prxmir_matched, int bph, int bpv, int & aln_size,
                     unsigned const & extend_direct, t_seq const & prx_seq, t_seq const & tag_seq, std::string const & org)
{
    bool align_on_prx = true;
    int num_match = 0;
    bool iso5p_canonic = false, iso3p_canonic = false;

    //debug_stream << tag_seq << "\n";
    //debug_stream << prx_seq << "\n";
    //debug_stream << bph << " : " << bpv << " : " << aln_size << "\n";
    if( extend_direct == EXTEND_LEFT || extend_direct == EXTEND_BOTH )
    {
        for(unsigned i = 1; (i <= bph) && (i <= bpv); ++i) // (i <= bpv) && ((bph - i) >= 0)
        {
            if(tag_seq[bpv - i] != prx_seq[bph - i]) {
                align_on_prx = false;
            }
            else {
                ++num_match;
            }
            //debug_stream << tag_seq[bpv - i] << " - " << prx_seq[bph - i] << "\t" << bpv - i << " - " << bph - i << "\n";
        }
        if(align_on_prx && ((bph - bpv ) >= 0)) {
            iso5p_canonic = true;
        }
    }
    if( extend_direct == EXTEND_RIGHT || extend_direct == EXTEND_BOTH )
    {
        for(unsigned i = 0; i < (tag_seq.size() - (bpv + aln_size) ) && ((bph + aln_size + i) < prx_seq.size()); ++i)
        {
            if(tag_seq[bpv + aln_size + i] != prx_seq[bph + aln_size + i]) {
                align_on_prx = false;
            }
            else {
                ++num_match;
            }
            //debug_stream << tag_seq[bpv + aln_size + i] << " - " << prx_seq[bph + aln_size + i] << "\t" << bpv + aln_size + i << " - " << bph + aln_size + i << "\n";
        }
        if(align_on_prx && ((bph + aln_size + num_match) <= prx_seq.size())) {
            iso3p_canonic = true;
        }
    }
    aln_size = aln_size + num_match; // here matches between tag and prxmir
    prxmir_matched.align_score = num_match;
    prxmir_matched.iso5p_canonic = iso5p_canonic;
    prxmir_matched.iso3p_canonic = iso3p_canonic;
}

// ----------------------------------------------------------------------------
// Function align_on_prxmir()
// ----------------------------------------------------------------------------

void align_on_prxmir(t_mir_matched_cell & mir_match, t_tag_cell & tag, t_mir_cell const & mir,
        unsigned const & extend_direct, int const & begin_diff, t_map_str_bool const & org_ids_m)
{
    if(extend_direct > 0) // if tag is not longer than mir we will not have any match with a best prxmir
    {
        int size_align = mir.seq.size();
        for(auto & org_prxmirs: mir.org_prx_m)
        {
            if(org_ids_m.count(ALL_ORGANISM) > 0 || org_ids_m.count(org_prxmirs.first) > 0) {
                for (unsigned i = 0; i < org_prxmirs.second.size(); ++i) {
                    t_prxmir_matched_cell prxmir_matched;
                    if (org_prxmirs.second[i].prx_seq.size() > 0) // if prxmir exist (maybe this check is useless)
                    {
                        //debug_stream << size_align << "\t" << org_prxmirs.first << "\n";
                        t_mir_info_cell mirna = t_mir_info_cell();
                        int start_gen;

                        if (mir.index == org_prxmirs.second[i].index5p) { mirna = org_prxmirs.second[i].mirna5p; }
                        else if (mir.index == org_prxmirs.second[i].index3p) { mirna = org_prxmirs.second[i].mirna3p; }

                        if (org_prxmirs.second[i].start_gen_pri < 0) {// if pri coordinates not present => miRBase
                            start_gen = org_prxmirs.second[i].start_gen_pre;
                        }
                        else {
                            start_gen = org_prxmirs.second[i].start_gen_pri;
                        }

                        unsigned int start_off_mirna = mirna.start_gen - start_gen;
                        unsigned int end_off_mirna = mirna.end_gen - start_gen;
                        unsigned int offset_in_prxmir;

                        if (mirna.strand == '+' || mirna.strand == '?' ||
                            mirna.mimat_prx.empty()) // if we are in miRGeneDb or in miRBase (with strand +)
                        {
                            offset_in_prxmir = start_off_mirna;
                        } else // if we are in miRBase (with strand -)
                        {
                            offset_in_prxmir = org_prxmirs.second[i].prx_seq.size() - end_off_mirna;
                        }
                        select_prxmir_ext(mir_match, prxmir_matched, offset_in_prxmir, begin_diff, size_align,
                                          extend_direct, org_prxmirs.second[i].prx_seq, tag.seq, org_prxmirs.first);
                        //debug_stream << org_prxmirs.second[i].prx_seq << "\n" << mirna.info << "\t" << mir.seq << "\t"
                        //             << offset_in_prxmir << "\n";
                        if (mir.index ==
                            org_prxmirs.second[i].index5p) { prxmir_matched.mir_info = &org_prxmirs.second[i].mirna5p; }
                        else if (mir.index ==
                                 org_prxmirs.second[i].index3p) { prxmir_matched.mir_info = &org_prxmirs.second[i].mirna3p; }
                        prxmir_matched.prxmir = &org_prxmirs.second[i];
                        prxmir_matched.org = org_prxmirs.first;
                        prxmir_matched.mir_info_id = mir.index;
                        prxmir_matched.prxmir_id = i;
                        //debug_stream << prxmir_matched.mir_info->index << "\t" << prxmir_matched.mir_info->info << "\tPTR\t" << prxmir_matched.mir_info << "\n";
                        mir_match.mir_prxmir[org_prxmirs.first].push_back(prxmir_matched);
                        //debug_stream << mir_match.mir_prxmir[org_prxmirs.first].back().mir_info->index << "\t" << mir_match.mir_prxmir[org_prxmirs.first].back().mir_info->info << "\tPTR\t" << mir_match.mir_prxmir[org_prxmirs.first].back().mir_info << "\t" << org_prxmirs.first << "\n";
                        //debug_stream << org_prxmirs.second[i].index << "\t" << org_prxmirs.second[i].info << "\t" << offset_in_prxmir << "\t" << org_prxmirs.second[i].prx_seq << "\t" << size_align << "\n";
                    }
                }
                std::sort(mir_match.mir_prxmir[org_prxmirs.first].begin(),
                          mir_match.mir_prxmir[org_prxmirs.first].end(),
                          compare_by_score<t_prxmir_matched_cell>);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function reduce_cigar()
// ----------------------------------------------------------------------------

void reduce_cigar(std::string & cigar_red, std::string const & cigar_exp)
{
    int count = 1;
    for(unsigned i = 1; i < cigar_exp.length(); ++i)
    {
        if(cigar_exp[i] == cigar_exp[i - 1]) {
            ++count;
        }
        else
        {
            if(count == 1)
            {
                cigar_red.push_back(cigar_exp[i - 1]);
            } else
            {
                cigar_red += std::to_string(count);
                cigar_red.push_back(cigar_exp[i - 1]);
            }
            count = 1;
        }
    }
    if(count == 1)
    {
        cigar_red.push_back(cigar_exp[cigar_exp.length() -1]);
    } else
    {
        cigar_red += std::to_string(count);
        cigar_red.push_back(cigar_exp[cigar_exp.length() -1]);
    }
}

// ----------------------------------------------------------------------------
// Function check_aligned_letters()
// ----------------------------------------------------------------------------

bool check_aligned_letters(char const & nt)
{
    if(nt == 'A' || nt == 'C' || nt == 'G' || nt == 'T' || nt == 'U' || nt == 'a' || nt == 'c' || nt == 'g'
       || nt == 'u' || nt == 't') {
        return false;
    }
    else {
        return true;
    }
}

// ----------------------------------------------------------------------------
// Function update_align_flags()
// ----------------------------------------------------------------------------

void update_align_flags(t_mir_matched_cell & mir_match, t_seq const & mir_seq, options const & options)
{
    unsigned m = 0;
    unsigned cent_site = 0, supp_site = 0, comp_site = 0;
    for(unsigned i = 0; i < mir_match.align_iupac_mir.size(); ++i)
    {
        if(check_aligned_letters(mir_match.align_iupac_mir[i])) {
            ++m;
        }
        if(i == (options.seed_end + 1) && !check_aligned_letters(mir_match.align_iupac_mir[i])) {//options.seed_end + 1 is the offset site
            mir_match.off_site = true;
        }
        if(i > START_CENTRAL_SITE && i <= END_CENTRAL_SITE && cent_site < SIZE_CENTRAL_SITE)
        {
            if(!check_aligned_letters(mir_match.align_iupac_mir[i - 1]) && !check_aligned_letters(mir_match.align_iupac_mir[i])) {
                ++cent_site;
            }
            else {
                cent_site = 0;
            }
        }
        if(i > START_SUPP_SITE && i <= END_SUPP_SITE && supp_site < SIZE_SUPP_SITE)
        {
            if(!check_aligned_letters(mir_match.align_iupac_mir[i - 1]) && !check_aligned_letters(mir_match.align_iupac_mir[i])) {
                ++supp_site;
            }
            else {
                supp_site = 0;
            }
        }
        if(i > START_COMP_SITE && i <= END_COMP_SITE && comp_site < SIZE_COMP_SITE)
        {
            if(!check_aligned_letters(mir_match.align_iupac_mir[i - 1]) && !check_aligned_letters(mir_match.align_iupac_mir[i])) {
                ++comp_site;
            }
            else {
                comp_site = 0;
            }
        }
        if(mir_match.align_iupac_mir[i] == 'R' && mir_seq[i].to_char() == 'A') // A2I -> A in miR and G in tag
        {
            if(i >= options.seed_start && i <= options.seed_end) {
                mir_match.a2i_seed = true;
            }
            else {
                mir_match.a2i_out_seed = true;
            }
        }
    }
    if(m == 1)
    {
        mir_match.iso_snp = true;
        mir_match.iso_msnp = false;
    } else if (m > 1)
    {
        mir_match.iso_snp = false;
        mir_match.iso_msnp = true;
    }
    if((mir_match.iso5p == 0) && (mir_match.iso3p == 0) && !mir_match.iso_snp && !mir_match.iso_msnp) {
        mir_match.mir_exact = true;
    }
    if(cent_site >= SIZE_CENTRAL_SITE) {
        mir_match.central_site = true;
    }
    if(supp_site >= SIZE_SUPP_SITE) {
        mir_match.suppl_site = true;
    }
    if(comp_site >= SIZE_COMP_SITE) {
        mir_match.compens_site = true;
    }
}

// ----------------------------------------------------------------------------
// Function compute_align_score_and_cigar()
// ----------------------------------------------------------------------------

void compute_align_score_and_cigar(t_mir_matched_v & mir_tag_out, t_tag_cell & tag, t_mir_cell const & mir,
        int const & bph, int const & bpv, int const & aln_size, int const & score_seed_align,
        t_map_str_bool const & org_ids_m, options const & options)
{
    t_mir_matched_cell mir_match;
    int min_size_seq = std::min(tag.seq.size(), mir.seq.size());
    mir_match.align_score = min_size_seq - ( std::max( bph, bpv ) * INDEL_PENALTY ) + ( score_seed_align * MISMATCH_ALIGN_PENALTY_SEED );
    mir_match.mismatch_in_seed = abs(score_seed_align);
    mir_match.align_begin_mir = bph;
    mir_match.align_begin_tag = bpv;
    mir_match.align_size = aln_size;
    mir_match.mir_tag_size_diff =  mir.seq.size() - tag.seq.size();
    std::string cigar;
    int begin_diff = bpv - bph;
    unsigned align_on_prxmir_code = 0;
    mir_match.iso5p = begin_diff;
    for(int j = 0; j < begin_diff; ++j) // 5p tag size > 5p mir size
    {
        if(j == 0)
        {
            align_on_prxmir_code = EXTEND_LEFT; // align on the prxmir list
        }
        mir_match.align_score = mir_match.align_score - INDEL_PENALTY;
        //debug_stream << "X\t - \t" << tag.seq[j] << "\t:\t" << mir_match.align_score << "\n";
        cigar.push_back('I');
        mir_match.align_iupac_tag.push_back(std::tolower(tag.seq[j].to_char()));
    }
    for(int i =  0; i < bph; ++i)
    {
        if((begin_diff + i) < 0)  // 5p tag size < 5p mir size
        {
            mir_match.align_score = mir_match.align_score - INDEL_PENALTY;
            //debug_stream << mir.seq[i] << "\t - \tX" << "\t:\t" << mir_match.align_score << "\n";
            cigar.push_back('D');
            mir_match.align_iupac_mir.push_back(std::tolower(mir.seq[i].to_char()));
        } else  // 5p tag size unaligned from 5p mir size
        {
            if(mir.seq[i] != tag.seq[i + begin_diff])
            {
                cigar.push_back(tag.seq[i+begin_diff].to_char());
                mir_match.align_score = mir_match.align_score - MISMATCH_UNALIGN_PENALTY;
            } else //TODO verify if is required a penalty even for the matching case
            {
                cigar.push_back('M');
            }
            //debug_stream << mir.seq[i] << "\t - \t" << tag.seq[i + begin_diff] << "\t:\t" << mir_match.align_score << "\n";
            mir_match.align_iupac_tag.push_back(nt_iupac_converter(mir.seq[i].to_char(), tag.seq[i+begin_diff].to_char()));
            mir_match.align_iupac_mir.push_back(nt_iupac_converter(mir.seq[i].to_char(), tag.seq[i+begin_diff].to_char()));
        }
    }
    for(int i = bph; i < aln_size; ++i) // tag aligned to mir
    {
        if(mir.seq[i] != tag.seq[i + begin_diff])
        {
            cigar.push_back(tag.seq[i+begin_diff].to_char());
            mir_match.align_score = mir_match.align_score - MISMATCH_ALIGN_PENALTY;
        } else //TODO verify if is required a penalty even for the matching case
        {
            cigar.push_back('M');
        }
        //debug_stream << mir.seq[i] << "\t - \t" << tag.seq[i + begin_diff] << "\t:\t" << mir_match.align_score << "\n";
        mir_match.align_iupac_tag.push_back(nt_iupac_converter(mir.seq[i].to_char(), tag.seq[i+begin_diff].to_char()));
        mir_match.align_iupac_mir.push_back(nt_iupac_converter(mir.seq[i].to_char(), tag.seq[i+begin_diff].to_char()));
    }
    for(int i = aln_size; i < mir.seq.size(); ++i)
    {
        if( (begin_diff + i) < tag.seq.size() ) // 3p tag size unaligned from 3p mir size
        {
            if(mir.seq[i] != tag.seq[i + begin_diff])
            {
                mir_match.align_score = mir_match.align_score - MISMATCH_UNALIGN_PENALTY;
                cigar.push_back(tag.seq[i+begin_diff].to_char());
            } else //TODO verify if is required a penalty even for the matching case
            {
                cigar.push_back('M');
            }
            //debug_stream << mir.seq[i] << "\t - \t" << tag.seq[i + begin_diff] << "\t:\t" << mir_match.align_score << "\n";
            mir_match.align_iupac_tag.push_back(nt_iupac_converter(mir.seq[i].to_char(), tag.seq[i+begin_diff].to_char()));
            mir_match.align_iupac_mir.push_back(nt_iupac_converter(mir.seq[i].to_char(), tag.seq[i+begin_diff].to_char()));
        }else // 3p tag size < 3p mir size
        {
            --mir_match.iso3p;
            mir_match.align_score = mir_match.align_score - INDEL_PENALTY;
            //debug_stream << mir.seq[i] << "\t-\tX" << "\t:\t" << mir_match.align_score << "\n";
            cigar.push_back('D');
            mir_match.align_iupac_mir.push_back(std::tolower(mir.seq[i].to_char()));
        }
    }
    for(int j = mir.seq.size() + begin_diff; j < tag.seq.size(); ++j) // 3p tag size > 3p mir size
    {
        if(j == ( mir.seq.size() + begin_diff ) ) // align on the prxmir list
        {
            if (align_on_prxmir_code == 0) {
                align_on_prxmir_code = EXTEND_RIGHT;
            }
            else {
                align_on_prxmir_code = EXTEND_BOTH;
            }
        }
        ++mir_match.iso3p;
        mir_match.align_score = mir_match.align_score - INDEL_PENALTY;
        //debug_stream << "X\t-\t" << tag.seq[j] << "\t:\t" << mir_match.align_score << "\n";
        cigar.push_back('I');
        mir_match.align_iupac_tag.push_back(std::tolower(tag.seq[j].to_char()));
    }
    if(mir_match.align_score >= options.min_align_score) // this will prevent to save and compute weak alignment
    {
        align_on_prxmir(mir_match, tag, mir, align_on_prxmir_code, begin_diff, org_ids_m);
        mir_match.mir = &mir; // collecting sequences
        reduce_cigar(mir_match.cigar, cigar);
/*
        debug_stream << cigar << "\n";
        debug_stream << mir_match.cigar << "\n";
        debug_stream << mir_match.align_iupac_mir << "\n";
        debug_stream << mir_match.align_iupac_tag << "\n";
        debug_stream << mir_match.align_score << "\n";

        if (mir_match.mir_prxmir.size() > 0)
        {
            for(unsigned w = 0; w < mir_match.mir_prxmir.at("hsa").size(); ++w)
            {
                debug_stream <<  mir_match.mir_prxmir.at("hsa").size() << "\n";
                debug_stream << mir_match.mir->seq << "\t" << mir_match.mir->index << "\n";
                debug_stream << tag.seq << "\t" << tag.index << "\n";
                debug_stream << "PAS_" << mir_match.mir_prxmir.at("hsa")[w].align_score << "\tIC5_" << mir_match.mir_prxmir.at("hsa")[w].iso5p_canonic
                          << "\tIC3_" << mir_match.mir_prxmir.at("hsa")[w].iso3p_canonic << "\n";
                debug_stream << "PTR\t" << mir_match.mir_prxmir.at("hsa")[w].mir_info << "\n";
                debug_stream << "\tMII_" << mir_match.mir_prxmir.at("hsa")[w].mir_info->index << "\n";
                debug_stream << "\tMRF_" << mir_match.mir_prxmir.at("hsa")[w].mir_info->ref << "\n";
                debug_stream << "\tMIN_" << mir_match.mir_prxmir.at("hsa")[w].mir_info->info[0] << "\n";
                debug_stream << "\tMRF_" << mir_match.mir_prxmir.at("hsa")[w].mir_info->ref << "\tMSG_" << mir_match.mir_prxmir.at("hsa")[w].mir_info->start_gen
                          << "\tMEG_" << mir_match.mir_prxmir.at("hsa")[w].mir_info->end_gen << "\tMMT_" << mir_match.mir_prxmir.at("hsa")[w].mir_info->mimat
                          << "\tMSR_" << mir_match.mir_prxmir.at("hsa")[w].mir_info->strand;
                debug_stream << "\tPII_" << mir_match.mir_prxmir.at("hsa")[w].prxmir->index << "\tPIN_" << mir_match.mir_prxmir.at("hsa")[w].prxmir->info
                          << "\tPRF_" << mir_match.mir_prxmir.at("hsa")[w].prxmir->ref << "\n";
            }
        }
*/
        update_align_flags(mir_match, mir.seq, options);
        mir_tag_out.push_back(mir_match); // add an alignment between mir and tag
/*
        if (mir_match.mir_prxmir.size() > 0)
        {
            for(unsigned w = 0; w < mir_match.mir_prxmir.at("hsa").size(); ++w)
            {
                debug_stream << "PAS_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].align_score
                             << "\tIC5_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].iso5p_canonic
                             << "\tIC3_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].iso3p_canonic << "\n";
                debug_stream << "\tMII_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].mir_info->index
                             << "\tMIN_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].mir_info->info
                             << "\tMRF_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].mir_info->ref
                             << "\tMSG_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].mir_info->start_gen
                             << "\tMEG_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].mir_info->end_gen
                             << "\tMMT_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].mir_info->mimat
                             << "\tMSR_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].mir_info->strand
                             << "\tPII_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].prxmir->index
                             << "\tPIN_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].prxmir->info
                             << "\tPRF_" << mir_tag_out.back().mir_prxmir.at("hsa")[w].prxmir->ref << "\n";
            }
        }
*/
    }
}

// ----------------------------------------------------------------------------
// Function ungapped_extension()
// ----------------------------------------------------------------------------

void ungapped_extension(t_mir_matched_v & mir_tag_out, t_tag_cell & tag, t_mir_cell const & mir, int const & begin_pos_tag_seed,
                       int const & end_pos_tag_seed, int score_seed_align, t_map_str_bool const & org_ids_m, options const & options)
{
    int bph{options.seed_start}, bpv{begin_pos_tag_seed}, aln_size{end_pos_tag_seed - begin_pos_tag_seed - 1}; // this -1 is for compensating the exclusion of last index in pattern match
    seed_ext_both(bph, bpv, aln_size, mir.seq, tag.seq); // with this extension we start from the seed found and align until up to 2 mismatches (one on the left and one on the right)
    seed_ext_both(bph, bpv, aln_size, mir.seq, tag.seq); // again up to 2 mismatches one per side
    if(aln_size > options.min_size_aln_stp1)
    {
        //debug_stream << mir.index << "\t" << mir.seq << "\n" << tag.index << "\t" << tag.seq << "\n";
        //debug_stream << begin_pos_tag_seed << " : " << end_pos_tag_seed << " = " << end_pos_tag_seed - begin_pos_tag_seed << "\n";
        //debug_stream << bph << " : " << bph + aln_size << " = " << aln_size << "\t" << mir.seq[bph] << " : " << mir.seq[bph + aln_size] << "\n";
        //debug_stream << bpv << " : " << bpv + aln_size << " = " << aln_size << "\t"  << tag.seq[bpv] << " : " << tag.seq[bpv + aln_size] << "\n";
        // Continue to align up to the maximum number of mismatches allowed
        for(int i = 0; i < options.mismatches_out_seed - 2; ++i)  // maximum number of allowed mismatches
        {
            seed_ext_rigth(bph, bpv, aln_size, mir.seq, tag.seq);
            //debug_stream << bph << " : " << bph + aln_size << " = " << aln_size << "\t" << mir.seq[bph] << " : " << mir.seq[bph + aln_size] << "\n";
            //debug_stream << bpv << " : " << bpv + aln_size << " = " << aln_size << "\t" << tag.seq[bpv] << " : " << tag.seq[bpv + aln_size] << "\n";
        }
        if(bpv + aln_size == tag.seq.size())
        {
            --aln_size;
            //debug_stream << bph << " : " << bph + aln_size << " = " << aln_size << "\t" << mir.seq[bph] << " : " << mir.seq[bph + aln_size] << "\n";
            //debug_stream << bpv << " : " << bpv + aln_size << " = " << aln_size << "\t" << tag.seq[bpv] << " : " << tag.seq[bpv + aln_size] << "\n";
        }
        if((aln_size >= options.min_size_tag) && (bpv <= options.max_start_pos_tag)  // alignments shorter than tag size are discarded
            && (tag.seq.size() <= (mir.seq.size() + (options.max_start_pos_tag * 2) ) ) ) // tag longer than mirna sequence + options.max_start_pos_tag * 2 are discarded
        {
            compute_align_score_and_cigar(mir_tag_out, tag, mir, bph, bpv, aln_size, score_seed_align, org_ids_m, options);
        }
    }
}

// ----------------------------------------------------------------------------
// Function scan_mir()
// ----------------------------------------------------------------------------

void scan_mir(t_mir_matched_v & mir_tag_out, t_tag_cell & tag, t_seed_v & seeds, t_map_str_bool const & org_ids_m, options const & options)
{
    seqan::IupacString tag_seqan2;
    from_rna15_to_iupac(tag_seqan2, tag.seq); // convert tag sequence (std::vector<rna15>) to IupacString in order to use seqan2 find algorithm
    t_prefix seed_zone_in_tag(tag_seqan2, options.max_start_pos_tag + options.seed_end - options.seed_start); // so we are safe that begin of the seed-tag alignment is lower than the user imposed threshold
    seqan::Finder<t_prefix> finder_tag4seed(seed_zone_in_tag);
    for(unsigned j = 0; j < seeds.size(); ++j)
    {
        seqan::goBegin(finder_tag4seed); // move Finder to the beginning of the text
        seqan::clear(finder_tag4seed); // reset Finder
        seqan::IupacString seed_seqan2;
        from_rna15_to_iupac(seed_seqan2, seeds[j].seq); // convert tag sequence (std::vector<rna15>) to IupacString in order to use seqan2 find algorithm
        seqan::Pattern<seqan::IupacString, seqan::Myers<>> pattern(seed_seqan2, -options.mismatches_in_seed);
        while (seqan::find(finder_tag4seed, pattern, -options.mismatches_in_seed))
        {
            while (seqan::findBegin(finder_tag4seed, pattern, seqan::getScore(pattern)))
            {
                if(seqan::length(seqan::infix(finder_tag4seed)) == (options.seed_end - options.seed_start + 1) ) // With this condition we select seeds that have only mismatches, if we want indels we should remove this condition
                {
                    //debug_stream << "SEED: {" << options.seed_start << ',' << options.seed_end << "]\t" << seeds[j].seq;
                    //debug_stream << "\tTAG: [" << seqan::beginPosition(finder_tag4seed) << ',' << seqan::endPosition(finder_tag4seed)
                    //             << ")\t" << seqan::infix(finder_tag4seed) << "\t" << seqan::length(seqan::infix(finder_tag4seed))
                    //             << "\t" << seqan::getScore(pattern) << "\t(tag: " << tag_seqan2 << ")\n";
                    for(unsigned k = 0; k < seeds[j].mir_seqs.size(); ++k)
                    {
                        /*bool mirna_ok = false;
                        if(org_ids_m.count(ALL_ORGANISM) <= 0) {
                            for (t_pair_str_bool const &org_bool : org_ids_m) {
                                if (!seeds[j].mir_seqs[k].org_prx_m[org_bool.first].empty()) {
                                    mirna_ok = true;
                                    break;
                                }
                            }
                        } else
                            mirna_ok = true;
                        if(mirna_ok) {*/
                            ungapped_extension(mir_tag_out, tag, seeds[j].mir_seqs[k],
                                               seqan::beginPosition(finder_tag4seed),
                                               seqan::endPosition(finder_tag4seed), seqan::getScore(pattern), org_ids_m,
                                               options);
                        //}
                    }
                }
            }
        }
    }
    std::sort(mir_tag_out.begin(), mir_tag_out.end(), compare_by_score<t_mir_matched_cell>);
}

#endif //ISOMIR_SEA_CORE_H
