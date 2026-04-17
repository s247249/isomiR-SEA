/*
 * =============================================================================
 * input.h
 *
 * Input file acquisition and processing.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_INPUT_H
#define ISOMIR_SEA_INPUT_H

#include "read_write_gff.h"

using namespace seqan3;

struct rna15_traits : sequence_file_input_default_traits_dna
{
    using sequence_alphabet = rna15; // instead of dna5
};

std::ostream& operator<<(std::ostream& output, const t_seq& seq)
{
    for(rna15 nucl : seq) { output << nucl.to_char(); }
    return output;
}

// ----------------------------------------------------------------------------
// Function extract_elements()
// ----------------------------------------------------------------------------

void extract_elements(t_map_str_bool & msb, std::string const & s, options const & options)
{
    t_vect_str tmpvs;
    //std::string delims = "!,#,$,%,&,',(,),*,+,,,-,.,/,:,;,<,=,>,?,[,\\,],^,_,";
    //boost::split(tmpvs, s, boost::is_any_of(delims));
    boost::split(tmpvs, s, is_punct);
    for(unsigned i = 0; i < tmpvs.size(); ++i)
    {
        std::transform(tmpvs[i].begin(), tmpvs[i].end(), tmpvs[i].begin(), [](unsigned char c){ return std::tolower(c); });
        msb[tmpvs[i]] = true;
    }
}

// ----------------------------------------------------------------------------
// Function filter_organism()
// ----------------------------------------------------------------------------

void filter_organism(t_vect_str & ids, std::vector<t_seq> & seqs, t_qual & quals, t_map_str_bool const & org_ids_m)
{
    t_vect_str ids_new;
    std::vector<t_seq> seqs_new;
    t_qual quals_new;

    for(unsigned i = 0; i < ids.size(); ++i)
    {
        std::string id_words, first_word;

        id_words = ids[i].substr(0, ids[i].find(' '));
        first_word = id_words.substr(0, id_words.find('-'));
        std::transform(first_word.begin(), first_word.end(), first_word.begin(), [](unsigned char c){ return std::tolower(c); }); // org_code
        if((org_ids_m.count(ALL_ORGANISM) > 0) || org_ids_m.count(first_word) > 0)
        {
            ids_new.push_back(ids[i]);
            seqs_new.push_back(seqs[i]);
            quals_new.push_back(quals[i]);
        }/*
        else // if you want to be more efficient (but does not preserve the order of the lines of the input file)
        {
            ids[i] = ids.back();
            ids.pop_back();
            seqs[i] = seqs.back();
            seqs.pop_back();
            quals[i] = quals.back();
            quals.pop_back();
            --i;
        }*/
    }
    ids = ids_new;
    seqs = seqs_new;
    quals = quals_new;
}

// ----------------------------------------------------------------------------
// Function acquire_from_fasta()
// ----------------------------------------------------------------------------

void acquire_from_fasta(t_vect_str & ids, std::vector<t_seq> & seqs, t_qual & quals, std::string const & filename)
{
    sequence_file_input<rna15_traits> seq_file_in{filename};
    for (auto & [seq, id, qual] : seq_file_in)
    {
        ids.push_back(id.c_str());
        seqs.push_back(seq);
        quals.push_back(qual);
    }
}

// ----------------------------------------------------------------------------
// Function acquire_from_fasta()
// ----------------------------------------------------------------------------

void acquire_from_fasta(t_vect_str & ids, std::vector<t_seq> & seqs, t_qual & quals, t_map_str_bool const & org_ids_m,
        std::string const & filename)
{
    acquire_from_fasta(ids, seqs, quals, filename);
    //filter_organism(ids, seqs, quals, org_ids_m);
}

// ----------------------------------------------------------------------------
// Function identify_premir_specs()
// ----------------------------------------------------------------------------

void identify_prxmir_specs(t_prx_mirna_cell & prxmir_cell, std::string & org_code, std::string const & id,
                           unsigned const & prxmir_index)
{
    t_vect_str id_words;
    std::string centr_word;

    boost::split(id_words, id, boost::is_any_of(" "));
    org_code = id_words[0].substr(0, id_words[0].find('-'));
    std::transform(org_code.begin(), org_code.end(), org_code.begin(), [](unsigned char c){ return std::tolower(c); }); // org_code
    prxmir_cell.index = prxmir_index;
    prxmir_cell.info = id; // full string reported in the reference DB
    centr_word = id_words[0].substr(0, id_words[0].find('_')); // this will remove from miRGeneDB the pri/pre label
    prxmir_cell.mir_id = centr_word.substr(ORGCODELENGTH + 1, centr_word.length()); // premiRNA id
    if(id_words.size() > 1) {
        prxmir_cell.mimat = id_words[1];  // mimat
    }
}

// ----------------------------------------------------------------------------
// Function fill_prxmir_vector()
// ----------------------------------------------------------------------------

void fill_prxmir_vector(t_org_prxmir_m & ref_prxmir_db, t_vect_str & ids, std::vector<t_seq> & seqs, t_qual & quals,
                        unsigned & prxmir_index)
{
    unsigned i = 0;
    for(; i < ids.size(); ++i) // i is the prxmir_index
    {
        t_prx_mirna_cell prxmir_cell;
        std::string org_code;

        prxmir_cell.prx_seq = seqs[i];
        identify_prxmir_specs(prxmir_cell, org_code, ids[i], i + prxmir_index);
        std::transform(prxmir_cell.mir_id.begin(), prxmir_cell.mir_id.end(), prxmir_cell.mir_id.begin(), [](unsigned char c){ return std::tolower(c); });
        if(prxmir_cell.mimat.empty()) {// miRGeneDb
            ref_prxmir_db[org_code][prxmir_cell.mir_id] = prxmir_cell;
        }
        else {// miRBase
            ref_prxmir_db[org_code][prxmir_cell.mimat] = prxmir_cell;
        }
    }
    prxmir_index = i + prxmir_index;
}

// ----------------------------------------------------------------------------
// Function acquire_from_gff()
// ----------------------------------------------------------------------------

void acquire_from_gff(std::vector<gff_record> & gff_records_mir, std::vector<gff_record> & gff_records_prxmir,
                      std::string const & filename, options const & options)
{
    std::ifstream in(filename);
    //auto stream_view = view::istreambuf(in);
    auto stream_view = std::ranges::subrange(
        std::istreambuf_iterator<char>(in),
        std::istreambuf_iterator<char>()
    );
    gff_record record;
    auto it = stream_view.begin();
    auto end = stream_view.end();
    while( it != end )
    {
        read_gff_record(record, it, end);
        if(record.type == options.type_mir_gff) // TODO verify this condition that can be week for some ref DBs
        {
            gff_records_mir.push_back(record);
        } else
        {
            gff_records_prxmir.push_back(record);
        }
        //debug_stream << record.ref << "\t" << record.source << "\t" << record.type << "\t" << record.begin_pos << "\t" << record.end_pos << "\t" << record.score << "\t" << record.strand << "\t" << record.phase << "\t";
        //for(int i=0; i<record.tag_names.size(); i++) { debug_stream << record.tag_names[i] << "=" << record.tag_values[i] << ";"; }
        //debug_stream << "\n";
    }
}

// ----------------------------------------------------------------------------
// Function identify_mir_specs()
// ----------------------------------------------------------------------------

void identify_mir_specs(t_mir_info_cell & mir_info_cell, std::string & org_code, std::string const & id, unsigned const & mir_index)
{
    t_vect_str id_words, last_word;
    std::string sdrow_di;

    boost::split(id_words, id, boost::is_any_of(" "));
    org_code = id_words[0].substr(0, id_words[0].find('-'));
    std::transform(org_code.begin(), org_code.end(), org_code.begin(), [](unsigned char c){ return std::tolower(c); }); // org_code
    mir_info_cell.info = id; // full string reported in the reference DB
    sdrow_di = id_words[0];
    std::reverse(sdrow_di.begin(), sdrow_di.end());
    boost::split(last_word, sdrow_di, boost::is_any_of("_"));
    if(last_word.size() < 2) // <2 means that there was not a "_" and so we are in miRBase // TODO Generalize this check (works with the current miRBase and MiRGeneDB) or see TODO in load_reference_db
    {
        last_word.clear();
        last_word.push_back(sdrow_di.substr(0, sdrow_di.find('-')));
    }
    bool tmp_flag = false;
    bool p53 = false;
    bool sm = false;
    unsigned it;
    unsigned of7 = 0;
    for(it = 0; it < last_word[0].length(); ++it)
    {
        if(last_word[0][it] == 'p')
        {
            tmp_flag = true;
            of7 += 2; //_p or -p are considered
        }
        if(last_word[0][it] == '5' && tmp_flag)
        {
            mir_info_cell.end5p3p = END5P; // flag 5p
            ++of7;
            p53 = true;
        }
        else if(last_word[0][it] == '3' && tmp_flag)
        {
            mir_info_cell.end5p3p = END3P; // flag 3p
            ++of7;
            p53 = true;
        }
        else if(last_word[0][it] == '*')
        {
            mir_info_cell.mature_star = STAR; // flag *
            ++of7;
            sm = true;
        }
    }
    mir_info_cell.index = mir_index;
    unsigned dim = 0;
    if(p53 || sm) {
        dim = id_words[0].length() - of7;
    }
    else {
        dim = id_words[0].length();
    }
    mir_info_cell.mir_id.clear();
    for(unsigned j = ORGCODELENGTH + 1; j < dim; ++j)
    {
        mir_info_cell.mir_id.push_back((id_words[0][j])); //miRNA id
    }
    if(id_words.size() > 1) {
        mir_info_cell.mimat = id_words[1];  //MIMAT
    }
};

// ----------------------------------------------------------------------------
// Function fill_mir_info_cell()
// ----------------------------------------------------------------------------

void fill_mir_info_cell(t_mir_info_cell & mirna, t_mir_info_cell const & mir_info_cell, gff_record & gff_record_mir,
                        std::string const & id, std::string const & tmp_id)
{
    mirna.info = mir_info_cell.info;
    mirna.mir_id = id;
    mirna.end5p3p = mir_info_cell.end5p3p;
    mirna.mature_star = mir_info_cell.mature_star;
    mirna.index = mir_info_cell.index;
    mirna.start_gen = gff_record_mir.begin_pos;
    mirna.end_gen = gff_record_mir.end_pos;
    mirna.ref = gff_record_mir.ref;
    mirna.strand = gff_record_mir.strand;
    if(mirna.mimat.empty() && gff_record_mir.tag_values.size() > 1) // TODO Generalize this check (works with the current miRBase and MiRGeneDB) or see TODO in load_reference_db
    {
        mirna.mimat = gff_record_mir.tag_values[1];
        if(gff_record_mir.tag_values.size() > 2) {
            mirna.mimat_prx = tmp_id;
        }
    }
}

// ----------------------------------------------------------------------------
// Function record_other_prxmir_location()
// ----------------------------------------------------------------------------

void record_other_prxmir_location(std::string & id, t_prx_mirna_cell & prxmir)
{
    t_oth_loc location;
    location.name = id;
    location.start_gen_pre = prxmir.start_gen_pre;
    location.end_gen_pre = prxmir.end_gen_pre;
    location.ref = prxmir.ref;
    prxmir.other_locations.push_back(location);
}

// ----------------------------------------------------------------------------
// Function update_gff_prxmir_records()
// ----------------------------------------------------------------------------

void update_gff_prxmir_records(std::vector<gff_record> & gff_records_prxmir, t_map_str_bool const & org_ids_m,
                               t_org_prxmir_m & ref_prxmir_db)
{
    for(unsigned i = 0; i < gff_records_prxmir.size(); ++i)
    {
        t_mir_info_cell mir_info_cell;
        std::string org_code, tmp_id;
        unsigned mirgene_mirbase_f = 0; // 0 = mirgene | 1 = mirbase

        if(gff_records_prxmir[i].tag_values.size() > 2) // miRBase // TODO Generalize this check (works with the current miRBase and MiRGeneDB) or see TODO in load_reference_db
        {
            identify_mir_specs(mir_info_cell, org_code, gff_records_prxmir[i].tag_values[2], i);
            //debug_stream << org_code << "\t" << mir_info_cell.mir_id << "\t" << mir_info_cell.end5p3p << std::endl;
            mir_info_cell.mimat = gff_records_prxmir[i].tag_values[0];
            mir_info_cell.mimat = mir_info_cell.mimat.substr(0, mir_info_cell.mimat.find('_'));
            mirgene_mirbase_f = 1;
        }
        else // miRGeneDB
        {
            identify_mir_specs(mir_info_cell, org_code, gff_records_prxmir[i].tag_values[0], i);
            std::string id = mir_info_cell.mir_id.substr(0, mir_info_cell.mir_id.find('_'));
            mir_info_cell.mir_id = id;
            //std::cout << org_code << "\t" << mir_info_cell.mir_id << "\t" << mir_info_cell.end5p3p << std::endl;
        }
        //if((org_ids_m.count(ALL_ORGANISM) > 0) || org_ids_m.count(org_code) > 0)
        //{
            std::transform(mir_info_cell.mir_id.begin(), mir_info_cell.mir_id.end(), mir_info_cell.mir_id.begin(), [](unsigned char c){ return std::tolower(c); }); // org_code
            if(mirgene_mirbase_f == 1) {
                tmp_id = mir_info_cell.mimat;
            }
            else {
                tmp_id = mir_info_cell.mir_id;
            }

            if(ref_prxmir_db[org_code][tmp_id].start_gen_pre >= 0)
            {
                // This prxmir occur more times in the gff file, every time with different genome location.
                // We store all the info of the last occurrence and only the genome location for the rest of them.
                // TODO Maybe review and implement a data structure that manages and saves all these different occurrences
                record_other_prxmir_location(tmp_id, ref_prxmir_db[org_code][tmp_id]);
            }
            ref_prxmir_db[org_code][tmp_id].start_gen_pre = gff_records_prxmir[i].begin_pos;
            ref_prxmir_db[org_code][tmp_id].end_gen_pre = gff_records_prxmir[i].end_pos;
            ref_prxmir_db[org_code][tmp_id].ref = gff_records_prxmir[i].ref;
            ref_prxmir_db[org_code][tmp_id].strand = gff_records_prxmir[i].strand;
            if (empty(ref_prxmir_db[org_code][tmp_id].mimat) && gff_records_prxmir[i].tag_values.size() > 1) {// TODO Generalize this check (works with the current miRBase and MiRGeneDB) or see TODO in load_reference_db
                ref_prxmir_db[org_code][tmp_id].mimat = gff_records_prxmir[i].tag_values[1];
            }
            if(gff_records_prxmir[i].tag_values.size() <= 2) // miRGeneDB // TODO Generalize this check (works with the current miRBase and MiRGeneDB) or see TODO in load_reference_db
            {
                ref_prxmir_db[org_code][tmp_id].start_gen_pri = gff_records_prxmir[i].begin_pos - 30; // 30-Flank
                ref_prxmir_db[org_code][tmp_id].end_gen_pri = gff_records_prxmir[i].end_pos + 30; // 30-Flank
            }
        //}
    }
}

// ----------------------------------------------------------------------------
// Function update_gff_mir_records()
// ----------------------------------------------------------------------------

void update_gff_mir_records(std::vector<gff_record> & gff_records_mir, t_map_str_bool const & org_ids_m,
                            t_org_prxmir_m & ref_prxmir_db)
{
    for(unsigned i = 0; i < gff_records_mir.size(); ++i)
    {
        t_mir_info_cell mir_info_cell;
        std::string org_code, tmp_id;
        unsigned mirgene_mirbase_f = 0;

        if(gff_records_mir[i].tag_values.size() > 2) // miRBase // TODO Generalize this check (works with the current miRBase and MiRGeneDB) or see TODO in load_reference_db
        {
            identify_mir_specs(mir_info_cell, org_code, gff_records_mir[i].tag_values[2], i);
            //debug_stream << org_code << "\t" << mir_info_cell.mir_id << "\t" << mir_info_cell.end5p3p << std::endl;
            mir_info_cell.mimat = gff_records_mir[i].tag_values[0];
            mir_info_cell.mimat = mir_info_cell.mimat.substr(0, mir_info_cell.mimat.find('_'));
            mir_info_cell.mimat_prx = gff_records_mir[i].tag_values[3];
            mirgene_mirbase_f = 1;
        }
        else // miRGeneDB
        {
            identify_mir_specs(mir_info_cell, org_code, gff_records_mir[i].tag_values[0], i);
            std::string id = mir_info_cell.mir_id.substr(0, mir_info_cell.mir_id.find('_'));
            mir_info_cell.mir_id = id;
            //std::cout << org_code << "\t" << mir_info_cell.mir_id << "\t" << mir_info_cell.end5p3p << std::endl;
        }
        //if((org_ids_m.count(ALL_ORGANISM) > 0) || org_ids_m.count(org_code) > 0)
        //{
            std::string id = mir_info_cell.mir_id; // To preserve letter case
            std::transform(mir_info_cell.mir_id.begin(), mir_info_cell.mir_id.end(), mir_info_cell.mir_id.begin(), [](unsigned char c){ return std::tolower(c); }); // org_code
            if(mirgene_mirbase_f == 1) {
                tmp_id = mir_info_cell.mimat_prx;
            }
            else {
                tmp_id = mir_info_cell.mir_id;
            }

            if(mir_info_cell.end5p3p == END3P)
            {
                fill_mir_info_cell(ref_prxmir_db[org_code][tmp_id].mirna3p, mir_info_cell, gff_records_mir[i], id, tmp_id);
            }
            else if(mir_info_cell.end5p3p == END5P)
            {
                fill_mir_info_cell(ref_prxmir_db[org_code][tmp_id].mirna5p, mir_info_cell, gff_records_mir[i], id, tmp_id);
            }
            else // ENDUNK
            {
                // Temporary save (or overwrite each time you encounter the same ID) the ENDUNK miRNA
                fill_mir_info_cell(ref_prxmir_db[org_code][tmp_id].mirna_tmp, mir_info_cell, gff_records_mir[i], id, tmp_id);
                ++ref_prxmir_db[org_code][tmp_id].mirna_tmp.count_overwrite;
                if(ref_prxmir_db[org_code][tmp_id].mirna_tmp.count_overwrite >= ref_prxmir_db[org_code][tmp_id].other_locations.size()+1)
                {
                    // This is the last occurrence of the ENDUNK miRNA
                    // We can understand if it is 5p or 3p looking if it is in the first half or second of the prxmir
                    int length_prxmir = ref_prxmir_db[org_code][tmp_id].end_gen_pre - ref_prxmir_db[org_code][tmp_id].start_gen_pre;
                    if(gff_records_mir[i].begin_pos - ref_prxmir_db[org_code][tmp_id].start_gen_pre < length_prxmir/2) // miRNA begin in the first half of the prxmir
                    {
                        if(gff_records_mir[i].strand == '+' || gff_records_mir[i].strand == '?' || mir_info_cell.mimat_prx.empty()) // if we are in miRGeneDb or in miRBase (with strand +)
                        {
                            mir_info_cell.info += "-5p";
                            mir_info_cell.end5p3p = END5P;
                            fill_mir_info_cell(ref_prxmir_db[org_code][tmp_id].mirna5p, mir_info_cell, gff_records_mir[i], id, tmp_id);
                        }
                        else
                        {
                            mir_info_cell.info += "-3p";
                            mir_info_cell.end5p3p = END3P;
                            fill_mir_info_cell(ref_prxmir_db[org_code][tmp_id].mirna3p, mir_info_cell, gff_records_mir[i], id, tmp_id);
                        }
                    }
                    else // mirna begin in the second half of the prxmir
                    {
                        if(gff_records_mir[i].strand == '+' || gff_records_mir[i].strand == '?' || mir_info_cell.mimat_prx.empty()) // if we are in miRGeneDb or in miRBase (with strand +)
                        {
                            mir_info_cell.info += "-3p";
                            mir_info_cell.end5p3p = END3P;
                            fill_mir_info_cell(ref_prxmir_db[org_code][tmp_id].mirna3p, mir_info_cell, gff_records_mir[i], id, tmp_id);
                        }
                        else
                        {
                            mir_info_cell.info += "-5p";
                            mir_info_cell.end5p3p = END5P;
                            fill_mir_info_cell(ref_prxmir_db[org_code][tmp_id].mirna5p, mir_info_cell, gff_records_mir[i], id, tmp_id);
                        }
                    }
                }
            }
        //}
    }
}

// ----------------------------------------------------------------------------
// Function update_gff_in_prxmir()
// ----------------------------------------------------------------------------

void update_gff_in_prxmir_vector(std::vector<gff_record> & gff_records_mir, std::vector<gff_record> & gff_records_prxmir,
                                 t_map_str_bool const & org_ids_m, t_org_prxmir_m & ref_prxmir_db)
{
    update_gff_prxmir_records(gff_records_prxmir, org_ids_m, ref_prxmir_db);
    update_gff_mir_records(gff_records_mir, org_ids_m, ref_prxmir_db);
}

// ----------------------------------------------------------------------------
// Function extract_sub_string()
// ----------------------------------------------------------------------------

void extract_sub_string(t_seq & out_sub_str, t_seq const & str, unsigned const & start, unsigned const & end,
        std::string const & mirna_info, options const & options)
{
    if(start < end && str.size() > end)
    {
        out_sub_str.resize(end - start + 1);
        for(unsigned i = start; i <= end; ++i)
        {
            out_sub_str[i - start] = str[i];
        }
    } else
    {
        std::string message = std::string("WARNING: error in substring extraction, ") +
                                          "probably the .gff file contains inconsistent coordinates with respect to " +
                                          "primary-sequences file, check: " + mirna_info;
        _VVV(options, message);
        out_sub_str = mirna_N;
    }
}

// ----------------------------------------------------------------------------
// Function find_element_or_create()
// ----------------------------------------------------------------------------

template <typename TElements, typename TElement>
void find_element_or_create(int & pos, TElements & xs, TElement & x)
{
    if(xs.size() > 0)
    {
        // Check if this seed is already in the vector
        for(unsigned j = 0; j < xs.size(); ++j)
        {
            if (xs[j].seq == x.seq) {
                pos = j;
            }
        }
    }
    if(pos < 0)
    {
        // Push a new element
        pos = xs.size();
        xs.push_back(x);
    }
}

// ----------------------------------------------------------------------------
// Function max_mir_size_update()
// ----------------------------------------------------------------------------

void max_mir_size_update(t_seq const & seq)
{
    int size = seq.size();
    if (size > MAX_MIR_SIZE) {
        MAX_MIR_SIZE = size;
    }
}

// ----------------------------------------------------------------------------
// Function fsv_using_mir_info_cell(), "fsv" means "fill seed vector"
// ----------------------------------------------------------------------------

void fsv_using_mir_info_cell(t_seed_v & seeds, t_prx_mirna_cell & prxmir, unsigned start_gen_pre, t_mir_info_cell & mirna,
                             unsigned & seq_index, options const & options)
{
    t_seed_cell seed_cell;
    t_mir_cell mir_cell;
    int seed_pos = -1, mir_pos = -1, record_position = -1;
    unsigned old_size;
    std::string info_words, org_code;
    unsigned int start_off_mirna = mirna.start_gen - start_gen_pre;
    unsigned int end_off_mirna = mirna.end_gen - start_gen_pre;

    // Extract mirna sequence
    if(mirna.strand == '+' || mirna.strand == '?' || mirna.mimat_prx.empty()) {// if we are in miRGeneDb or in miRBase (with strand +)
        extract_sub_string(mir_cell.seq, prxmir.prx_seq, start_off_mirna, end_off_mirna - 1, mirna.info, options);
    }
    else {// if we are in miRBase (with strand -)
        extract_sub_string(mir_cell.seq, prxmir.prx_seq, prxmir.prx_seq.size() - end_off_mirna,
                           prxmir.prx_seq.size() - start_off_mirna - 1, mirna.info, options);
    }
    //debug_stream << prxmir.prx_seq << "\n" << mirna.info << " " << mir_cell.seq << "\n";
    max_mir_size_update(mir_cell.seq);

    // Extract seed
    extract_sub_string(seed_cell.seq, mir_cell.seq, options.seed_start, options.seed_end, mirna.info, options);
    //debug_stream << seed_cell.seq << "\n";

    find_element_or_create(seed_pos, seeds, seed_cell);
    seeds[seed_pos].index = seed_pos;
    old_size = seeds[seed_pos].mir_seqs.size();
    find_element_or_create(mir_pos, seeds[seed_pos].mir_seqs, mir_cell);
    if(mir_pos == old_size) { // if not present, assign new index
        seeds[seed_pos].mir_seqs[mir_pos].index = seq_index;
        record_position = seq_index;
    } else { // if present retrieve old index
        record_position = seeds[seed_pos].mir_seqs[mir_pos].index;
    }
    if(mirna.end5p3p == END5P) {
        prxmir.index5p = record_position;
    } else {
        prxmir.index3p = record_position;
    }
    seq_index = seq_index + seeds[seed_pos].mir_seqs.size() - old_size;
    //debug_stream << seeds[seed_pos].mir_seqs[mir_pos].seq << " : " << seeds[seed_pos].mir_seqs[mir_pos].index << "\n";
    info_words = mirna.info.substr(0, mirna.info.find(' '));
    org_code = info_words.substr(0, info_words.find('-'));
    std::transform(org_code.begin(), org_code.end(), org_code.begin(), [](unsigned char c){ return std::tolower(c); }); // org_code
    seeds[seed_pos].mir_seqs[mir_pos].org_prx_m[org_code].push_back(prxmir);
}

// ----------------------------------------------------------------------------
// Function fill_seed_vector()
// ----------------------------------------------------------------------------

void fill_seed_vector(t_seed_v & seeds, t_org_prxmir_m & ref_prxmir_db, unsigned seq_index, options const & options)
{
    for(auto & org: ref_prxmir_db)
    {
        for (auto & prxmir: org.second)
        {
            int start_gen;
            if(prxmir.second.start_gen_pri < 0) {// if pri coordinates not present => miRBase
                start_gen = prxmir.second.start_gen_pre;
            }
            else {
                start_gen = prxmir.second.start_gen_pri;
            }
            if(prxmir.second.mirna5p.index >= 0) {
                fsv_using_mir_info_cell(seeds, prxmir.second, start_gen, prxmir.second.mirna5p, seq_index, options);
            }
            if(prxmir.second.mirna3p.index >= 0) {
                fsv_using_mir_info_cell(seeds, prxmir.second, start_gen, prxmir.second.mirna3p, seq_index, options);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function fill_mir_info_cell_from_mature()
// ----------------------------------------------------------------------------

void fill_mir_info_cell_from_mature(t_mir_info_cell & mirna, t_vect_str & id_words, unsigned const & END,
        unsigned & mirna_index, bool & sm)
{
    mirna.info = id_words[0];
    for(int k=2; k<id_words.size()-1; k++)
    {
        mirna.info.append(" ");
        mirna.info.append(id_words[k]);
    }
    mirna.mir_id = id_words[0];
    mirna.end5p3p = END;
    if(sm) {
        mirna.mature_star = STAR;
    }
    mirna.index = mirna_index;
    if(id_words.size() > 1)
    {
        mirna.mimat = id_words[1];
    }
}

// ----------------------------------------------------------------------------
// Function identify_prxmir_and_mir_specs()
// ----------------------------------------------------------------------------

void identify_prxmir_and_mir_specs(t_prx_mirna_cell & prxmir_cell, std::string & org_code, std::string const & id,
                                   t_seq const & seq, unsigned const & prxmir_index, unsigned & mirna_index)
{ // TODO Eventually generalize this function (works with the current miRBase and MiRGeneDB) or see TODO in load_reference_db
    t_vect_str id_words, last_word;
    std::string centr_word, sdrow_di;

    boost::split(id_words, id, boost::is_any_of(" "));
    org_code = id_words[0].substr(0, id_words[0].find('-'));
    std::transform(org_code.begin(), org_code.end(), org_code.begin(), [](unsigned char c){ return std::tolower(c); }); // org_code
    prxmir_cell.index = prxmir_index;
    sdrow_di = id_words[0];
    std::reverse(sdrow_di.begin(), sdrow_di.end());
    boost::split(last_word, sdrow_di, boost::is_any_of("_"));
    if(last_word.size() < 2) // <2 means that there was not a "_" and so we are in miRBase
    {
        last_word.clear();
        last_word.push_back(sdrow_di.substr(0, sdrow_di.find('-')));
    }
    bool tmp_flag = false;
    unsigned it;
    unsigned of7 = 0;
    bool p5 = false;
    bool p3 = false;
    bool sm = false;
    for(it = 0; it < last_word[0].length(); ++it)
    {
        if(last_word[0][it] == 'p')
        {
            tmp_flag = true;
            of7 += 2; //_p or -p are considered
        }
        if(last_word[0][it] == '5' && tmp_flag)
        {
            p5 = true; // flag 5p
            ++of7;
        }
        else if(last_word[0][it] == '3' && tmp_flag)
        {
            p3 = true; // flag 3p
            ++of7;
        }
        else if(last_word[0][it] == '*')
        {
            sm = true; // flag *
            ++of7;
        }
    }
    unsigned dim = 0;
    if(p5 || p3 || sm) {
        dim = id_words[0].length() - of7;
    }
    else {
        dim = id_words[0].length();
    }
    prxmir_cell.mir_id.clear();
    for(unsigned j = ORGCODELENGTH + 1; j < dim; ++j)
    {
        prxmir_cell.mir_id.push_back((id_words[0][j])); //miRNA id
        prxmir_cell.info.push_back((id_words[0][j]));
    }
    prxmir_cell.info.append(" ");
    for(int k=2; k<id_words.size()-1; k++)
    {
        prxmir_cell.info.append(id_words[k]);
        prxmir_cell.info.append(" ");
    }

    if(p5 || (!p5 && !p3)) // mirna5p (if 5p/3p is not specified, we set it to 5p)
    {
        fill_mir_info_cell_from_mature(prxmir_cell.mirna5p, id_words, END5P, mirna_index, sm);
        prxmir_cell.tmp5p = seq; // temporary store 5p here
    }
    else if(p3) // mirna3p
    {
        fill_mir_info_cell_from_mature(prxmir_cell.mirna3p, id_words, END3P, mirna_index, sm);
        prxmir_cell.tmp3p = seq; // temporary store 3p here
    }
    ++mirna_index;
}

// ----------------------------------------------------------------------------
// Function check_gff_correctness()
// ----------------------------------------------------------------------------

void check_gff_correctness(t_prx_mirna_cell & prxmir, t_seq & seq)
{
    // Works only with miRGeneDB sequences
    auto it_found = std::search(prxmir.prx_seq.begin(), prxmir.prx_seq.end(), seq.begin(), seq.end());
    int start_off_mat = it_found - prxmir.prx_seq.begin();
    int end_off_mat = start_off_mat + seq.size();
    int start_gen;

    if(prxmir.start_gen_pri < 0) {// if pri coordinates not present => miRBase
        start_gen = prxmir.start_gen_pre;
    }
    else {
        start_gen = prxmir.start_gen_pri;
    }

    if(start_off_mat < prxmir.prx_seq.size()/2) // 5p
    {
        int start_off_mirna = prxmir.mirna5p.start_gen - start_gen;
        int end_off_mirna = prxmir.mirna5p.end_gen - start_gen;

        // If coordinates differs, update them
        if(start_off_mat != start_off_mirna || end_off_mat != end_off_mirna) {
            //debug_stream << prxmir.info << "\n" << start_off_mat << " " << end_off_mat << "\n" << start_off_mirna << " " << end_off_mirna << "\n";
            prxmir.mirna5p.start_gen = start_off_mat + start_gen;
            prxmir.mirna5p.end_gen = end_off_mat + start_gen;
        }
    } else // 3p
    {
        int start_off_mirna = prxmir.mirna3p.start_gen - start_gen;
        int end_off_mirna = prxmir.mirna3p.end_gen - start_gen;

        // If coordinates differs, update them
        if(start_off_mat != start_off_mirna || end_off_mat != end_off_mirna) {
            //debug_stream << prxmir.info << "\n" << start_off_mat << " " << end_off_mat << "\n" << start_off_mirna << " " << end_off_mirna << "\n";
            prxmir.mirna3p.start_gen = start_off_mat + start_gen;
            prxmir.mirna3p.end_gen = end_off_mat + start_gen;
        }
    }
}

// ----------------------------------------------------------------------------
// Function fill_prxmirs_and_mirs_from_mature()
// ----------------------------------------------------------------------------

void fill_prxmirs_and_mirs_from_mature(t_org_prxmir_m & ref_prxmir_db, t_vect_str & ids, std::vector<t_seq> & seqs,
        t_qual & quals, unsigned & prxmir_index, unsigned & mirna_index, options const & options, int & file_mature_mirbase)
{
    unsigned i = 0;
    for(; i < ids.size(); ++i)
    {
        t_prx_mirna_cell prxmir_cell;
        std::string org_code;

        identify_prxmir_and_mir_specs(prxmir_cell, org_code, ids[i], seqs[i], prxmir_index, mirna_index);
        std::transform(prxmir_cell.mir_id.begin(), prxmir_cell.mir_id.end(), prxmir_cell.mir_id.begin(), [](unsigned char c){ return std::tolower(c); });
        if(!options.in_file_primir.empty() && (!prxmir_cell.mirna5p.mimat.empty() || !prxmir_cell.mirna5p.mimat.empty())) {
            file_mature_mirbase = 1;
            return;
        }
        if(ref_prxmir_db[org_code][prxmir_cell.mir_id].index >= 0) // prxmir_cell already created
        {
            if(!ref_prxmir_db[org_code][prxmir_cell.mir_id].prx_seq.empty()) // prxmir_cell already created with primir file
            {
                check_gff_correctness(ref_prxmir_db[org_code][prxmir_cell.mir_id], seqs[i]);
            }
            else  // prxmir_cell already created with mature/star file (5p or 3p already present)
            {
                // TODO eventually handle the case in which for a single mirID there is both the real 5p/3p and the same sequence without 5p/3p label
                if(ref_prxmir_db[org_code][prxmir_cell.mir_id].tmp3p.empty())
                {
                    ref_prxmir_db[org_code][prxmir_cell.mir_id].mirna3p = prxmir_cell.mirna3p;
                    ref_prxmir_db[org_code][prxmir_cell.mir_id].tmp3p = prxmir_cell.tmp3p;
                }
                else if(ref_prxmir_db[org_code][prxmir_cell.mir_id].tmp5p.empty())
                {
                    ref_prxmir_db[org_code][prxmir_cell.mir_id].mirna5p = prxmir_cell.mirna5p;
                    ref_prxmir_db[org_code][prxmir_cell.mir_id].tmp5p = prxmir_cell.tmp5p;
                }
            }
        }
        else // prxmir_cell not present
        {
            ref_prxmir_db[org_code][prxmir_cell.mir_id] = prxmir_cell;
            ++prxmir_index;
        }
    }
}

// ----------------------------------------------------------------------------
// Function build_prxmirs()
// ----------------------------------------------------------------------------

void build_prxmirs(t_org_prxmir_m & ref_prxmir_db, options const & options)
{
    for(auto & org: ref_prxmir_db)
    {
        for (auto & prxmir: org.second)
        {
            if(prxmir.second.prx_seq.empty())
            {
                prxmir.second.start_gen_pri = 0;
                prxmir.second.start_gen_pre = FLANK_SIZE;

                if(prxmir.second.tmp5p.empty())
                {
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), flank_N.begin(), flank_N.end());
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), mirna_N.begin(), mirna_N.end());
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), loop_N.begin(), loop_N.end());
                }
                else
                {
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), flank_N.begin(), flank_N.end());
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), prxmir.second.tmp5p.begin(), prxmir.second.tmp5p.end());
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), loop_N.begin(), loop_N.end());

                    prxmir.second.mirna5p.start_gen = FLANK_SIZE;
                    prxmir.second.mirna5p.end_gen = FLANK_SIZE + prxmir.second.tmp5p.size();

                    prxmir.second.tmp5p.clear();
                }
                if(prxmir.second.tmp3p.empty())
                {
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), mirna_N.begin(), mirna_N.end());
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), flank_N.begin(), flank_N.end());
                }
                else
                {
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), prxmir.second.tmp3p.begin(), prxmir.second.tmp3p.end());
                    prxmir.second.prx_seq.insert(prxmir.second.prx_seq.end(), flank_N.begin(), flank_N.end());

                    prxmir.second.mirna3p.start_gen = prxmir.second.prx_seq.size() - prxmir.second.tmp3p.size() - FLANK_SIZE;
                    prxmir.second.mirna3p.end_gen = prxmir.second.prx_seq.size() - FLANK_SIZE;

                    prxmir.second.tmp3p.clear();
                }
                prxmir.second.end_gen_pre = prxmir.second.prx_seq.size() - FLANK_SIZE;
                prxmir.second.end_gen_pri = prxmir.second.prx_seq.size();
                //debug_stream << prxmir.second.mir_id << " : " << prxmir.second.info << "\n";
                //debug_stream << prxmir.second.prx_seq << "\n";
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function print_mirna_info_cell()
// ----------------------------------------------------------------------------

template <typename t_file>
void print_mirna_info_cell(t_file & out_stream, t_mir_info_cell const & mirna)
{
    out_stream << "\t" << mirna.index << " : " << mirna.mir_id << " : "
               << mirna.mature_star << " : " << mirna.mimat << " : "
               << mirna.end5p3p << " : " << mirna.info
               << " : " << mirna.start_gen << " : " << mirna.end_gen
               << " : " << mirna.strand << " : " << mirna.ref << " : " << mirna.mimat_prx << "\n";
}

// ----------------------------------------------------------------------------
// Function print_other_locations()
// ----------------------------------------------------------------------------

template <typename t_file>
void print_other_locations(t_file & out_stream, std::vector<t_oth_loc> const & other_locations)
{
    out_stream << "Other prxmir locations: ";
    for(t_oth_loc loc : other_locations)
    {
        out_stream << loc.name << " : " << loc.start_gen_pre << " : " << loc.end_gen_pre << " : " << loc.ref << " | ";
    }
    out_stream << "\n";
}

// ----------------------------------------------------------------------------
// Function print_prxmir()
// ----------------------------------------------------------------------------
template <typename t_file>
void print_prxmir(t_file & out_stream, t_org_prxmir_m & ref_prxmir_db)
{
    for(t_org_prxmir_p const & org: ref_prxmir_db)
    {
        for (t_prx_mir_p const & prxmir: org.second)
        {
            out_stream << org.first << " : " << prxmir.second.index << " : " << prxmir.second.mir_id << " : "
                       << prxmir.first << " : " << prxmir.second.prx_seq << " : " << prxmir.second.mimat
                       << " : " << prxmir.second.start_gen_pre << " : " << prxmir.second.end_gen_pre
                       << " : " << prxmir.second.start_gen_pri << " : " << prxmir.second.end_gen_pri
                       << " : " << prxmir.second.strand << " : " << prxmir.second.ref << "\n";
            if(prxmir.second.mirna5p.index >= 0) {
                print_mirna_info_cell(out_stream, prxmir.second.mirna5p);
            }
            if(prxmir.second.mirna3p.index >= 0) {
                print_mirna_info_cell(out_stream, prxmir.second.mirna3p);
            }
            if(prxmir.second.other_locations.size() >= 0) {
                print_other_locations(out_stream, prxmir.second.other_locations);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function print_seeds()
// ----------------------------------------------------------------------------

template <typename t_file>
void print_seeds(t_file & out_stream, t_seed_v & seeds)
{
    for(unsigned i = 0; i < seeds.size(); ++i)
    {
        out_stream << seeds[i].index << "\t" << seeds[i].seq << "\n";
        for(unsigned j = 0; j < seeds[i].mir_seqs.size(); ++j)
        {
            out_stream << seeds[i].mir_seqs[j].index << "\t" << seeds[i].mir_seqs[j].seq << "\n";
            for(t_org_prx_p const & prxmir: seeds[i].mir_seqs[j].org_prx_m)
            {
                for (unsigned w = 0; w < prxmir.second.size(); ++w)
                {
                    out_stream << prxmir.first << " : " << prxmir.second[w].index << " : " << prxmir.second[w].mir_id
                               << " : " << prxmir.first << " : " << prxmir.second[w].prx_seq << " : " << prxmir.second[w].mimat
                               << " : " << prxmir.second[w].start_gen_pre << " : " << prxmir.second[w].end_gen_pre
                               << " : " << prxmir.second[w].start_gen_pri << " : " << prxmir.second[w].end_gen_pri
                               << " : " << prxmir.second[w].strand << " : " << prxmir.second[w].ref << "\n";
                    if(prxmir.second[w].mirna5p.index >= 0) {
                        print_mirna_info_cell(out_stream, prxmir.second[w].mirna5p);
                    }
                    if(prxmir.second[w].mirna3p.index >= 0) {
                        print_mirna_info_cell(out_stream, prxmir.second[w].mirna3p);
                    }
                    if(prxmir.second[w].other_locations.size() >= 0) {
                        print_other_locations(out_stream, prxmir.second[w].other_locations);
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function print_seed_mir_prx()
// ----------------------------------------------------------------------------

template <typename t_file>
void print_seed_mir_prx(t_file & out_stream, t_seed_v const & seeds)
{
    out_stream << "seedId" << "\tseedSeq" << "\tmirSeqId" << "\tmirSeqSeq" << "\torganism"
               << "\tmirInId" << "\tmirInName" << "\tmirInMatStar" << "\tmirInMimat" << "\tmirInEnd5p3p"
               << "\tmirInInfo" << "\tmirInStartGen" << "\tmirInEndGen" << "\tmirInStrand" << "\tmirInRef"
               << "\tprxId" << "\tprxInfo" << "\tprxMimat" << "\tprxRef" << "\tprxSeq"
               << "\tprxPreStartGen" << "\tprxPreEndGen" << "\tprxPriStartGen" << "\tprxPriEndGen"
               << "\tprxStrand" << "\n";
    for(unsigned i = 0; i < seeds.size(); ++i)
    {
        for(unsigned j = 0; j < seeds[i].mir_seqs.size(); ++j)
        {
            t_mir_cell const & mir_cell = seeds[i].mir_seqs[j];
            for(t_org_prx_p const & prxmir: mir_cell.org_prx_m)
            {
                for (unsigned w = 0; w < prxmir.second.size(); ++w)
                {
                    // Seed
                    out_stream << seeds[i].index << "\t" << seeds[i].seq << "\t";
                    // MirSeq
                    out_stream << mir_cell.index << "\t" << mir_cell.seq << "\t";
                    // Organism
                    out_stream << prxmir.first << "\t";
                    // MirInfo
                    t_mir_info_cell mirna;
                    if(mir_cell.index == prxmir.second[w].index5p) { mirna = prxmir.second[w].mirna5p; }
                    else if(mir_cell.index == prxmir.second[w].index3p) { mirna = prxmir.second[w].mirna3p; }
                    out_stream << mirna.index << "\t" << mirna.mir_id << "\t" << mirna.mature_star
                              << "\t" << mirna.mimat << "\t" << mirna.end5p3p << "\t" << mirna.info
                              << "\t" << mirna.start_gen << "\t" << mirna.end_gen
                              << "\t" << mirna.strand << "\t" << mirna.ref;
                    // PrxMir
                    out_stream << "\t" << prxmir.second[w].mir_id << "\t" << prxmir.second[w].info << "\t"
                               << prxmir.second[w].mimat << "\t" << prxmir.second[w].ref << "\t"
                               << prxmir.second[w].prx_seq << "\t"
                               << prxmir.second[w].start_gen_pre << "\t" << prxmir.second[w].end_gen_pre << "\t"
                               << prxmir.second[w].start_gen_pri << "\t" << prxmir.second[w].end_gen_pri << "\t"
                               << prxmir.second[w].strand << "\n";
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function load_reference_db()
// ----------------------------------------------------------------------------

void load_reference_db(t_seed_v & seeds, t_org_prxmir_m & ref_prxmir_db, t_map_str_bool const & org_ids_m, options const & options)
{ // TODO eventually the user must be informed about the format to use (because currently decisions are made based on the current version of miRBase and mirGeneDb databases)
    t_vect_str ids;
    std::vector<t_seq> seqs;
    t_qual quals;
    std::vector<gff_record> gff_records_mir, gff_records_prxmir;
    unsigned prxmir_index = 0, mirna_index = 0, seq_index = 0;
    int file_mature_mirbase = 0;

    format_fasta::file_extensions.push_back("fas"); // because ".fas" extension is not included in seqan3::format_fasta{}
    if(!options.unform_mature_names)  // this if we can read seqs from standard references such as miRBase or miRGeneDB
    {
        if(!options.in_file_primir.empty()) // if there is the primir file, there will also be the gff one
        {
            _VV(options, "STEP: Read reference database of primir sequences and gff coordinates.");
            acquire_from_fasta(ids, seqs, quals, org_ids_m, options.in_file_primir);
            fill_prxmir_vector(ref_prxmir_db, ids, seqs, quals, prxmir_index);
            //print_prxmir(debug_stream, ref_prxmir_db);

            acquire_from_gff(gff_records_mir, gff_records_prxmir, options.in_file_gff, options);
            update_gff_in_prxmir_vector(gff_records_mir, gff_records_prxmir, org_ids_m, ref_prxmir_db);
            //print_prxmir(debug_stream, ref_prxmir_db);
        }
        if(!options.in_file_mature.empty())
        {
            ids.clear();
            std::for_each(seqs.begin(), seqs.end(), [](auto& seq) { seq.clear(); }); seqs.clear();
            std::for_each(quals.begin(), quals.end(), [](auto& qual) { qual.clear(); }); quals.clear();
            _VV(options, "STEP: Read reference database of mature sequences.");
            acquire_from_fasta(ids, seqs, quals, org_ids_m, options.in_file_mature);
            fill_prxmirs_and_mirs_from_mature(ref_prxmir_db, ids, seqs, quals, prxmir_index, mirna_index, options, file_mature_mirbase);
            if(!options.in_file_star.empty())
            {
                ids.clear();
                std::for_each(seqs.begin(), seqs.end(), [](auto& seq) { seq.clear(); }); seqs.clear();
                std::for_each(quals.begin(), quals.end(), [](auto& qual) { qual.clear(); }); quals.clear();
                _VV(options, "STEP: Read reference database of star sequences.");
                acquire_from_fasta(ids, seqs, quals, org_ids_m, options.in_file_star);
                fill_prxmirs_and_mirs_from_mature(ref_prxmir_db, ids, seqs, quals, prxmir_index, mirna_index, options, file_mature_mirbase);
            }
            if(!file_mature_mirbase)
                build_prxmirs(ref_prxmir_db, options);
        }
        fill_seed_vector(seeds, ref_prxmir_db, seq_index, options);
        //print_seeds(debug_stream, seeds);
    } else
    {
        std::cerr << "File unformatted: must be designed a function capable to create a virtual management system for the data structure " <<  "\n"; //TODO implement this functionality
    }
}

// ----------------------------------------------------------------------------
// Function load_tags()
// ----------------------------------------------------------------------------

void load_tags(t_tag & tag, options const & options)
{
    std::ifstream indata;
    std::string line, word;
    unsigned tag_index = 0;

    indata.open(options.in_file_tags);
    while(std::getline(indata,line))
    {
        t_vect_str list_of_words;
        t_tag_cell tag_cell;

        boost::split(list_of_words, line, is_space);
        tag_cell.seq = list_of_words[0] | seqan3::views::char_to<seqan3::rna15> | ranges::to<std::vector<seqan3::rna15>>;
        if(list_of_words.size() > 2)
        {
            tag_cell.qual = list_of_words[1];
            tag_cell.count = atoi(list_of_words[2].c_str());
        }
        else
        {
            tag_cell.count = atoi(list_of_words[1].c_str());
        }
        tag_cell.index = tag_index;
        ++tag_index;
        if (list_of_words[0].size() > options.min_size_tag && ( list_of_words[0].size() <= ( MAX_MIR_SIZE + (options.max_start_pos_tag * 2) ) ) )
        {
            //debug_stream << tag_cell.seq << "\t" << tag_cell.qual << "\t" << tag_cell.count << "\t" << "\n";
            tag.read_count += tag_cell.count;
            tag.tag_v.push_back(tag_cell);
        } else // tag shorter than min_size_tag OR longer than mirna sequences + options.max_start_pos_tag * 2 are discarded
        {
            tag.invalid_read_count += tag_cell.count;
            tag.discarded_tag_v.push_back(tag_cell);
        }
    }
    tag.tag_count = tag.tag_v.size();
    tag.invalid_tag_count = tag.discarded_tag_v.size();
}

// ----------------------------------------------------------------------------
// Function load_tags_fastx()
// ----------------------------------------------------------------------------

void load_tags_fastx(t_tag & tag, options const & options)
{
    sequence_file_input<rna15_traits> seq_file_in{options.in_file_tags};
    for (auto & [seq, id, qual] : seq_file_in)
    {
        t_tag_cell tag_cell;
        t_vect_str id_fields;

        tag_cell.seq = seq;
        boost::split(id_fields, id, boost::is_any_of("|"));
        for(auto & field : id_fields) {
            if(field.rfind("ID:", 0) == 0)
            {
                tag_cell.index = std::stoul(field.substr(3, std::string::npos));
            }
            if(field.rfind("CN:", 0) == 0)
            {
                tag_cell.count = std::stoul(field.substr(3, std::string::npos));
            }
        }
        for(auto & qual_sign : qual){
            tag_cell.qual += qual_sign.to_char();
        }
        if (seq.size() > options.min_size_tag && ( seq.size() <= ( MAX_MIR_SIZE + (options.max_start_pos_tag * 2) ) ) )
        {
            tag.read_count += tag_cell.count;
            tag.tag_v.push_back(tag_cell);
        } else // tag shorter than min_size_tag OR longer than mirna sequences + options.max_start_pos_tag * 2 are discarded
        {
            tag.invalid_read_count += tag_cell.count;
            tag.discarded_tag_v.push_back(tag_cell);
        }
    }
    tag.tag_count = tag.tag_v.size();
    tag.invalid_tag_count = tag.discarded_tag_v.size();
}

void load_mirtrace_fastx(t_tag & tag, options const & options)
{
    sequence_file_input<rna15_traits> seq_file_in{options.in_file_tags};
    for (auto & [seq, id, qual] : seq_file_in)
    {
        t_tag_cell tag_cell;
        t_vect_str id_main;
        t_vect_str id_fields;

        tag_cell.seq = seq;
        
        // id = seq_1_x1234567 rnatype:mirna
        // id_main = { "seq_1_x1234567", "rnatype:mirna" }
        boost::split(id_main, id, boost::is_any_of(" "));
        // id_fields = { "seq", "1", "x1234567" }
        boost::split(id_fields, id_main[0], boost::is_any_of("_"));
        for (auto & field : id_fields) {
            if (field.rfind("x", 0) == 0)
            {
                tag_cell.count = std::stoul(field.substr(1));
            } else if (std::all_of(field.begin(), field.end(), ::isdigit))
            {
                tag_cell.index = std::stoul(field);
            }
        }

        for(auto & qual_sign : qual){
            tag_cell.qual += qual_sign.to_char();
        }
        if (seq.size() > options.min_size_tag && ( seq.size() <= ( MAX_MIR_SIZE + (options.max_start_pos_tag * 2) ) ) )
        {
            tag.read_count += tag_cell.count;
            tag.tag_v.push_back(tag_cell);
        } else // tag shorter than min_size_tag OR longer than mirna sequences + options.max_start_pos_tag * 2 are discarded
        {
            tag.invalid_read_count += tag_cell.count;
            tag.discarded_tag_v.push_back(tag_cell);
        }
    }
    tag.tag_count = tag.tag_v.size();
    tag.invalid_tag_count = tag.discarded_tag_v.size();
}


// ----------------------------------------------------------------------------
// Function prune_datastructure()
// ----------------------------------------------------------------------------

void prune_datastructure(t_seed_v & seeds, t_map_str_bool & org_ids_m)
{
    for(t_seed_cell & seed_cell : seeds) {
        seed_cell.mir_seqs.erase(std::remove_if(seed_cell.mir_seqs.begin(), seed_cell.mir_seqs.end(),
                [org_ids_m](t_mir_cell mir_cell){
                    bool mirna_ok = true;
                    if(org_ids_m.count(ALL_ORGANISM) <= 0) {
                        for (t_pair_str_bool const & org_bool : org_ids_m) {
                            if (!mir_cell.org_prx_m[org_bool.first].empty()) {
                                mirna_ok = false;
                                break;
                            }
                        }
                    } else
                        mirna_ok = false;
                    return mirna_ok;
        }), seed_cell.mir_seqs.end());
    }

    seeds.erase(std::remove_if(seeds.begin(), seeds.end(),
            [](t_seed_cell seed_cell){
                if(seed_cell.mir_seqs.size()<=0) { return true; }
                else { return false; }
    }), seeds.end());
}

#endif //ISOMIR_SEA_INPUT_H
