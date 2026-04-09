/*
 * =============================================================================
 * output.h
 *
 * Functions to generate output files.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_OUTPUT_H
#define ISOMIR_SEA_OUTPUT_H

namespace fs = std::experimental::filesystem;
typedef std::experimental::filesystem::path path;

// ============================================================================
// Functions
// ============================================================================

std::ostream& operator<<(std::ostream &out, const t_variant &v)
{
    std::visit([&out](auto&& arg) { out << arg; }, v);
    return out;
}

// ----------------------------------------------------------------------------
// Function bool2char()
// ----------------------------------------------------------------------------

char bool2char(bool const & b)
{
    if(b) { return 'T'; }
    else { return 'F'; }
}

// ----------------------------------------------------------------------------
// Function check_out_path()
// ----------------------------------------------------------------------------

bool check_out_path(options & options)
{
    if(options.out_file_label.length() == 0) // if common name empty
    {
        options.out_file_label = path(options.in_file_tags).filename().replace_extension("").generic_string(); //Use mirDB name as common part for the output file
        //_V(options, options.out_file_label);
    }
    if(options.path_out_files.length() > 0) // if output path assigned
    {
        if(path(options.path_out_files).empty())
        {
            std::cerr << "Invalid path of output file. Give a valid path or leave it empty for generating out in the tag file folder" << std::endl;
            return true;
        } else if (!is_directory(path(options.path_out_files)))
        {
            create_directory(path(options.path_out_files));
        } else if (path(options.path_out_files).is_relative())
        {
            options.path_out_files = path(options.path_out_files).relative_path().generic_string(); // if path exist, the absolute path is retrived and saved insead of relative path
        }
        //_V(options, options.path_out_files);
    } else
    {
        options.path_out_files = path(options.in_file_tags).parent_path().generic_string(); // if not assigned, path of the tag is used
        //_V(options, options.path_out_files);
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function print_option()
// ----------------------------------------------------------------------------

void print_options(std::ofstream & out_file, options const & options)
{
    out_file << "# Options of IsomiR-SEA run " << std::endl;
    out_file << "verbose = " << options.verbose << std::endl;
    out_file << "in_file_tags = " << options.in_file_tags << std::endl;
    out_file << "in_file_mature = " << options.in_file_mature << std::endl;
    out_file << "in_file_star = " << options.in_file_star << std::endl;
    out_file << "in_file_gff = " << options.in_file_gff << std::endl;
    out_file << "type_mir_gff = " << options.type_mir_gff << std::endl;
    //out_file << "in_file_bed = " << options.in_file_bed << std::endl;
    //out_file << "in_file_premir = " << options.in_file_premir << std::endl;
    out_file << "path_load_serialized = " << options.path_load_serialized << std::endl;
    out_file << "path_store_serialized = " << options.path_store_serialized << std::endl;
    out_file << "in_file_primir = " << options.in_file_primir << std::endl;
    out_file << "unform_mature_names = " << options.unform_mature_names << std::endl;
    out_file << "dna_or_rna = " << options.dna_or_rna << std::endl;
    out_file << "multi_sample_tab = " << options.multi_sample_tab << std::endl;
    out_file << "out_file_label = " << options.out_file_label << std::endl;
    out_file << "out_formats = " << options.out_formats << std::endl;
    out_file << "path_out_files = " << options.path_out_files << std::endl;
    out_file << "path_discarded_tags = " << options.path_discarded_tags << std::endl;
    out_file << "tmp_dir = " << options.tmp_dir << std::endl;

    //out_file << "dpallign_on_premirs = " << options.dp_allign_on_premirs << std::endl;
    //out_file << "affine_linear_dgs = " << options.affine_linear_dgs << std::endl;
    //out_file << "gap_open = " << options.gap_open << std::endl;
    //out_file << "gap_extend = " << options.gap_extend << std::endl;
    //out_file << "match_score = " << options.match_score << std::endl;
    //out_file << "mismatch_score = " << options.mismatch_score << std::endl;
    //out_file << "global_local = " << options.global_local << std::endl;
    //out_file << "un_top = " << options.un_top << std::endl;
    //out_file << "un_left = " << options.un_left << std::endl;
    //out_file << "un_right = " << options.un_right << std::endl;
    //out_file << "un_down = " << options.un_down << std::endl;

    out_file << "specie_codes = " << options.specie_codes << std::endl;
    out_file << "min_size_aln_stp1 = " << options.min_size_aln_stp1 << std::endl;
    out_file << "min_align_score = " << options.min_align_score << std::endl;
    out_file << "seed_start = " << options.seed_start << std::endl;
    out_file << "seed_end = " << options.seed_end << std::endl;
    out_file << "max_start_pos_tag = " << options.max_start_pos_tag << std::endl;
    out_file << "min_size_tag = " << options.min_size_tag << std::endl;
    out_file << "mismatches_out_seed = " << options.mismatches_out_seed << std::endl;
    out_file << "mismatches_in_seed = " << options.mismatches_in_seed << std::endl;

    out_file << "mir_db_name = " << options.mir_db_name << std::endl;
    out_file << "mir_db_version = " << options.mir_db_version << std::endl;
    out_file << "mir_db_web = " << options.mir_db_web << std::endl;

    out_file << "time_limit = " << options.time_limit << std::endl;
    out_file << "threads = " << options.threads << std::endl;
}

// ----------------------------------------------------------------------------
// Function fill_tag_mir_match_fields()
// ftmm = Fields Tag Mir (match) Map; mtoc = Mir Tag Out Cell; mtoc_mm = Mir Tag Out Cell MirMatch.
// ----------------------------------------------------------------------------

void fill_tag_mir_match_fields(t_map_str_any & ftmm, t_tag_mir_matched_v const & mtoc, t_mir_matched_cell const & mtoc_mm)
{
    // Tag
    ftmm[tmf[0].first] = std::make_pair(tmf[0].second, mtoc.tag->index);
    ftmm[tmf[1].first] = std::make_pair(tmf[1].second, mtoc.tag->seq);
    ftmm[tmf[2].first] = std::make_pair(tmf[2].second, mtoc.tag->qual);
    ftmm[tmf[3].first] = std::make_pair(tmf[3].second, mtoc.tag->count);
    // Mir
    ftmm[tmf[7].first] = std::make_pair(tmf[7].second, mtoc_mm.mir->index);
    ftmm[tmf[8].first] = std::make_pair(tmf[8].second, mtoc_mm.mir->seq);
    // Align
    ftmm[tmf[9].first] = std::make_pair(tmf[9].second, mtoc_mm.align_score);
    ftmm[tmf[10].first] = std::make_pair(tmf[10].second, mtoc_mm.align_iupac_tag);
    ftmm[tmf[11].first] = std::make_pair(tmf[11].second, mtoc_mm.align_iupac_mir);
    ftmm[tmf[12].first] = std::make_pair(tmf[12].second, mtoc_mm.cigar);
    ftmm[tmf[13].first] = std::make_pair(tmf[13].second, mtoc_mm.align_size);
    ftmm[tmf[14].first] = std::make_pair(tmf[14].second, mtoc_mm.mir_tag_size_diff);
    // Isomir
    ftmm[tmf[15].first] = std::make_pair(tmf[15].second, bool2char(mtoc_mm.mir_exact));
    ftmm[tmf[16].first] = std::make_pair(tmf[16].second, mtoc_mm.iso5p);
    ftmm[tmf[17].first] = std::make_pair(tmf[17].second, bool2char(mtoc_mm.iso_msnp));
    ftmm[tmf[18].first] = std::make_pair(tmf[18].second, bool2char(mtoc_mm.iso_snp));
    ftmm[tmf[19].first] = std::make_pair(tmf[19].second, mtoc_mm.iso3p);
    // Interaction Sites
    ftmm[tmf[20].first] = std::make_pair(tmf[20].second, bool2char(mtoc_mm.mismatch_in_seed));
    ftmm[tmf[21].first] = std::make_pair(tmf[21].second, bool2char(mtoc_mm.off_site));
    ftmm[tmf[22].first] = std::make_pair(tmf[22].second, bool2char(mtoc_mm.suppl_site));
    ftmm[tmf[23].first] = std::make_pair(tmf[23].second, bool2char(mtoc_mm.compens_site));
    ftmm[tmf[24].first] = std::make_pair(tmf[24].second, bool2char(mtoc_mm.central_site));
    // A2I
    ftmm[tmf[25].first] = std::make_pair(tmf[26].second, bool2char(mtoc_mm.a2i_seed));
    ftmm[tmf[26].first] = std::make_pair(tmf[26].second, bool2char(mtoc_mm.a2i_out_seed));
}

// ----------------------------------------------------------------------------
// Function fill_tag_mir_info_match_fields()
// ----------------------------------------------------------------------------

void fill_tag_mir_info_match_fields(t_map_str_any & ftmm, t_mir_info_cell const & mir_info, int const & mir_tag_size_diff,
        int const & tag_size)
{
    ftmm[tmf[29].first] = std::make_pair(tmf[29].second, mir_info.index);
    ftmm[tmf[30].first] = std::make_pair(tmf[30].second, mir_info.info);
    ftmm[tmf[31].first] = std::make_pair(tmf[31].second, mir_info.ref);
    ftmm[tmf[32].first] = std::make_pair(tmf[32].second, mir_info.start_gen);
    ftmm[tmf[33].first] = std::make_pair(tmf[33].second, mir_info.end_gen);
    ftmm[tmf[34].first] = std::make_pair(tmf[34].second, mir_info.mimat);
    ftmm[tmf[35].first] = std::make_pair(tmf[35].second, mir_info.strand);
    ftmm[tmf[4].first] = std::make_pair(tmf[4].second, mir_info.start_gen + mir_tag_size_diff); // tag coordinate start
    ftmm[tmf[5].first] = std::make_pair(tmf[5].second, mir_info.start_gen + mir_tag_size_diff + tag_size); // tag coordinate end
}

// ----------------------------------------------------------------------------
// Function fill_tag_prxmir_info_fields()
// ----------------------------------------------------------------------------

void fill_tag_prxmir_info_fields(t_map_str_any & ftmm, t_prx_mirna_cell const & prxmir)
{
    ftmm[tmf[36].first] = std::make_pair(tmf[36].second, prxmir.index);
    ftmm[tmf[37].first] = std::make_pair(tmf[37].second, prxmir.info);
    ftmm[tmf[38].first] = std::make_pair(tmf[38].second, prxmir.ref);
    ftmm[tmf[39].first] = std::make_pair(tmf[39].second, prxmir.prx_seq);
    if(prxmir.start_gen_pri < 0)
    {
        ftmm[tmf[40].first] = std::make_pair(tmf[40].second, prxmir.start_gen_pre);
        ftmm[tmf[41].first] = std::make_pair(tmf[41].second, prxmir.end_gen_pre);
    }
    else
    {
        ftmm[tmf[40].first] = std::make_pair(tmf[40].second, prxmir.start_gen_pri);
        ftmm[tmf[41].first] = std::make_pair(tmf[41].second, prxmir.end_gen_pri);
    }
    ftmm[tmf[42].first] = std::make_pair(tmf[42].second, prxmir.mimat);
    ftmm[tmf[43].first] = std::make_pair(tmf[43].second, prxmir.strand);
}

class my_size //: public boost::static_visitor<int>
{
public:
    int operator()(int i) const {
        return 1;
    }
    int operator()(unsigned i) const {
        return 1;
    }
    int operator()(char c) const {
        return 1;
    }
    int operator()(const std::string & str) const {
        return str.length();
    }
    int operator()(const t_seq & seq) const {
        return seq.size();
    }
};

class my_string //: public boost::static_visitor<std::string>
{
public:
    std::string operator()(int i) const
    {
        return std::to_string(i);
    }
    std::string operator()(unsigned i) const
    {
        return std::to_string(i);
    }
    std::string operator()(char c) const
    {
        std::string str(1, c);
        return str;
    }
    std::string operator()(const std::string & str) const
    {
        return str;
    }
    std::string operator()(const t_seq & seq) const
    {
        std::string str;
        for(seqan3::rna15 nucl : seq) { str += nucl.to_char(); }
        return str;
    }
};

// ----------------------------------------------------------------------------
// Function print_tag_mir_tab()
// ----------------------------------------------------------------------------

template <typename t_paths>
void print_tag_mir_tab(t_paths & tm_files, t_map_str_any const & ftmm, std::ifstream & multisample_table,
        std::string & line, int & line_n, options const & options)
{
    // Print alignment lines
    if(options.multi_sample_tab.empty()) // only one .tab output file
    {
        //for(unsigned i = 0; i < num_align; ++i) {
            for(unsigned j = 0; j < tmf.size(); ++j) {
                if(ftmm.count(tmf[j].first) > 0 && std::visit(my_size(), ftmm.at(tmf[j].first).second) > 0) {
                    tm_files[0] << ftmm.at(tmf[j].first).second << "\t";
                }
                else {
                    tm_files[0] << "?" << "\t"; //TODO verify if ? is appropriate for downstream analysis
                }
            }
            tm_files[0] << "\n";
        //}
    }
    else
    {
        t_vect_str columns;
        //for(unsigned i = 0; i < num_align; ++i) {
            unsigned tag_id = std::stoul(std::visit(my_string(), ftmm.at(tmf[0].first).second));
            columns.clear();

            while(tag_id != line_n) // search the line in the tabular file correspondant to tag_id
            {
                line.clear();
                std::getline(multisample_table, line);
                ++line_n;
            }
            boost::split(columns, line, boost::is_any_of("\t"));
            for(unsigned k=0; k<tm_files.size(); ++k) { // each columns is the # of reads of the tag in a specific sample
                int reads_in_k_file = std::stoi(columns[k + 1]); // k_file correspond to k+1_column
                if(reads_in_k_file > 0)
                {
                    for(unsigned j = 0; j < tmf.size(); ++j)
                    {
                        if(j!=3) // retrieve tag count (#reads) in a sample from the tabular file and not from tag-input-file
                        {
                            if(ftmm.count(tmf[j].first) > 0 && std::visit(my_size(), ftmm.at(tmf[j].first).second) > 0) {
                                tm_files[k] << ftmm.at(tmf[j].first).second << "\t";
                            }
                            else {
                                tm_files[k] << "?" << "\t"; //TODO verify if ? is appropriate for downstream analysis
                            }
                        }
                        else
                        {
                            tm_files[k] << reads_in_k_file << "\t";
                        }
                    }
                    tm_files[k] << "\n";
                }
            }
        //}
    }
}

// ----------------------------------------------------------------------------
// Function copy_var()
// ----------------------------------------------------------------------------

template <typename t_var_out, typename t_var_in>
void copy_var(t_var_out & var_out, t_var_in const & var_in)
{
    if(std::visit(my_size(), var_in) > 0) { var_out = std::get<t_var_out>(var_in); }
}

// ----------------------------------------------------------------------------
// Function append_gff_tag()
// ----------------------------------------------------------------------------

template <typename t_gff_tag, typename t_str, typename t_var_in>
void append_gff_tag(t_gff_tag & tag_names, t_gff_tag & tag_values, t_str const & var_in_id, t_var_in const & var_in_val)
{
    if(var_in_id.size() > 0)
    {
        tag_names.push_back(var_in_id);
        if(std::visit(my_size(), var_in_val) > 0)
        {
            tag_values.push_back(std::visit(my_string(), var_in_val));
        }
    }
}

// ----------------------------------------------------------------------------
// Function append_gff_tag2()
// ----------------------------------------------------------------------------

template <typename t_gff_tag, typename t_tmf>
void append_gff_tag2(t_gff_tag & tag_names, t_gff_tag & tag_values, t_map_str_any const & ftmm, t_tmf const & tmf)
{
    if(ftmm.count(tmf.first) > 0 && std::visit(my_size(), ftmm.at(tmf.first).second) > 0)
    {
        tag_names.push_back(ftmm.at(tmf.first).first);
        tag_values.push_back(std::visit(my_string(), ftmm.at(tmf.first).second));
    } else
    {
        tag_names.push_back(tmf.second);
        tag_values.push_back("?");
    }
}

// ----------------------------------------------------------------------------
// Function convert2gff()
// ----------------------------------------------------------------------------

void convert2gff(gff_record & record, t_map_str_any const & ftmm, bool const & first_mir, int & reads_in_k_file, options const & options)
{
    char tmp_c = '?';
    if(ftmm.count(tmf[31].first) > 0) {
        copy_var(record.ref, ftmm.at(tmf[31].first).second);
    }
    else if(ftmm.count(tmf[38].first) > 0) {
        copy_var(record.ref, ftmm.at(tmf[38].first).second);
    }
    record.source = options.mir_db_name; // TODO decide if keep as it is or add also the DB version
    if(ftmm.count(tmf[15].first) > 0)
    {
        copy_var(tmp_c, ftmm.at(tmf[15].first).second);
        if(tmp_c == 'T') {
            record.type = "ref_miRNA";
        }
        else if (tmp_c == 'F') {
            record.type = "isomiR";
        }
        else {
            record.type = "unk";
        }
    }
    else {
        record.type = "prxmiR";
    }
    int tmp_I;
    copy_var(tmp_I, ftmm.at(tmf[4].first).second);
    record.begin_pos = tmp_I;
    copy_var(tmp_I, ftmm.at(tmf[5].first).second);
    record.end_pos = tmp_I;
    copy_var(tmp_I, ftmm.at(tmf[9].first).second);
    record.score = tmp_I;
    if(ftmm.count(tmf[35].first) > 0) {
        copy_var(record.strand, ftmm.at(tmf[35].first).second);
    }
    else if(ftmm.count(tmf[43].first) > 0) {
        copy_var(record.strand, ftmm.at(tmf[43].first).second);
    }

    append_gff_tag2(record.tag_names, record.tag_values, ftmm, tmf[0]);
    append_gff_tag2(record.tag_names, record.tag_values, ftmm, tmf[1]);
    if(options.multi_sample_tab.empty()) {
        append_gff_tag2(record.tag_names, record.tag_values, ftmm, tmf[3]);
    } else { // retrieve tag count (#reads) in a sample from the tabular file and not from tag-input-file
        std::string tmp_name = "TC";
        t_variant tmp_variant = reads_in_k_file;
        append_gff_tag(record.tag_names, record.tag_values, tmp_name, tmp_variant);
    }
    append_gff_tag2(record.tag_names, record.tag_values, ftmm, tmf[37]);
    append_gff_tag2(record.tag_names, record.tag_values, ftmm, tmf[12]);
    if(ftmm.count(tmf[8].first) > 0)
    {
        append_gff_tag2(record.tag_names, record.tag_values, ftmm, tmf[30]);
        std::string iso, inter, tmp_name;
        t_variant tmp_variant;
        tmp_name = "ISO";
        iso = std::visit(my_string(), ftmm.at(tmf[15].first).second);
        if (ftmm.count(tmf[46].first) > 0) {
            iso += std::visit(my_string(), ftmm.at(tmf[46].first).second);
        }
        else {
            iso += "F";
        }
        iso = iso + std::visit(my_string(), ftmm.at(tmf[16].first).second)
              + std::visit(my_string(), ftmm.at(tmf[17].first).second)
              + std::visit(my_string(), ftmm.at(tmf[18].first).second)
              + std::visit(my_string(), ftmm.at(tmf[19].first).second);
        if (ftmm.count(tmf[47].first) > 0) {
            iso += std::visit(my_string(), ftmm.at(tmf[47].first).second);
        }
        else {
            iso += "F";
        }
        tmp_variant = iso;
        append_gff_tag(record.tag_names, record.tag_values, tmp_name, tmp_variant);
        tmp_name = "INT";
        inter = std::visit(my_string(), ftmm.at(tmf[20].first).second)
                + std::visit(my_string(), ftmm.at(tmf[21].first).second)
                + std::visit(my_string(), ftmm.at(tmf[22].first).second)
                + std::visit(my_string(), ftmm.at(tmf[23].first).second)
                + std::visit(my_string(), ftmm.at(tmf[24].first).second);
        tmp_variant = inter;
        append_gff_tag(record.tag_names, record.tag_values, tmp_name, tmp_variant);
        tmp_name = "NotPass";
        if(first_mir)
        {
            if ( (ftmm.count(tmf[27].first) > 0) && (std::visit(my_string(), ftmm.at(tmf[27].first).second) == "1") ) // this is not necessary, since following checks consider these
            {
                tmp_name = "Pass";
            }
            else if (ftmm.count(tmf[28].first) > 0)
            {
                if ((std::visit(my_string(), ftmm.at(tmf[28].first).second) != "0"
                     && std::visit(my_string(), ftmm.at(tmf[28].first).second) != "?")) {
                    tmp_name = "Pass";
                }
            }
        }
        tmp_variant = tmp_name;
        tmp_name = "FILTER";
        append_gff_tag(record.tag_names, record.tag_values, tmp_name, tmp_variant);
    } else
    {/*
        // Entering in an unused zone of code;
        std::string tmp_name;
        t_variant tmp_variant;
        //if (std::visit(my_string(), ftmm.at(tmf[44].first).second) == "1") //TODO improve this check for sequences that share the same scores
        //  tmp_name = "Pass";
        //else
        //  tmp_name = "NotPass";
        tmp_name = "Unknown";
        tmp_variant = tmp_name;
        tmp_name = "FILTER";
        append_gff_tag(record.tag_names, record.tag_values, tmp_name, tmp_variant);
    */}
    // TODO if necessary eventually implement a fully compliant mirTop notation
}

// ----------------------------------------------------------------------------
// Function print_tag_mir_gff()
// ----------------------------------------------------------------------------

template <typename t_paths>
void print_tag_mir_gff(t_paths & tm_files, t_map_str_any const & ftmm, std::ifstream & multisample_table,
        std::string & line, bool & first_mir, options const & options)
{
    t_vect_str columns;

    //for(unsigned i = 0; i < num_align; ++i)
    //{
        if(options.multi_sample_tab.empty()) // only one .gff output file
        {
            gff_record record;
            int tmp = -1;
            convert2gff(record, ftmm, first_mir, tmp, options);
            write_gff_record(tm_files[0], record);
        }
        else
        {
            unsigned tag_id = std::stoul(std::visit(my_string(), ftmm.at(tmf[0].first).second));
            columns.clear();

            boost::split(columns, line, boost::is_any_of("\t"));
            for(unsigned k=0; k<tm_files.size(); ++k) { // each columns is the # of reads of the tag in a specific sample
                gff_record record;
                int reads_in_k_file = std::stoi(columns[k + 1]); // k_file correspond to k+1_column
                if (reads_in_k_file > 0) {
                    convert2gff(record, ftmm, first_mir, reads_in_k_file, options);
                    write_gff_record(tm_files[k], record);
                }
            }
        }
    //}
}

// ----------------------------------------------------------------------------
// Function single_multi_discarded_old()
// ----------------------------------------------------------------------------
/*
void single_multi_discarded_old(t_map_str_any & ftmm, std::string & TI_old, bool & discard_all, bool & multi_mapped, bool & same_mir_diff_prx)
{
    int code = 3; // single-mapped: 1; multi-mapped: 2; discarded: 3;

    if(TI_old == "" || TI_old != (std::visit(my_string(), ftmm.at(tmf[0].first).second))) { // new tag
        discard_all = false;
        multi_mapped = false;
        if ((ftmm.count(tmf[27].first) > 0) && (std::visit(my_string(), ftmm.at(tmf[27].first).second) == "1")) {
            code = 1;
        }
        else if(ftmm.count(tmf[28].first) > 0)
        {
            if(std::visit(my_string(), ftmm.at(tmf[28].first).second) != "0"
                 && std::visit(my_string(), ftmm.at(tmf[28].first).second) != "?") {
                code = 1;
                if(same_mir_diff_prx) { code = 2; }
            }
            else if(std::visit(my_string(), ftmm.at(tmf[28].first).second) == "0") {
                if(ftmm.count(tmf[48].first) > 0 && std::visit(my_string(), ftmm.at(tmf[48].first).second) != "0"
                                                    && std::visit(my_string(), ftmm.at(tmf[48].first).second) != "?") {
                    code = 1;
                   discard_all = true;
                }
                else {
                    code = 2;
                    multi_mapped = true;
                }
            }
        }
    }
    else { // same tag
        if(ftmm.count(tmf[28].first) > 0 && !discard_all) {
            if (std::visit(my_string(), ftmm.at(tmf[28].first).second) != "0"
                 && std::visit(my_string(), ftmm.at(tmf[28].first).second) != "?") {
                if(multi_mapped) { code = 2; multi_mapped = false; } // last alignment with same score of multi-mapped is still multi-mapped
                else if(same_mir_diff_prx) { code = 2; }
                else { code = 3; }
                discard_all = true;
            }
            else if (std::visit(my_string(), ftmm.at(tmf[28].first).second) == "0") {
                code = 2;
            }
        }
        else {
            if(multi_mapped) { code = 2; multi_mapped = false; } // last alignment with same score of multi-mapped is still multi-mapped
            else { code = 3; }
        }
    }
    ftmm[tmf[54].first] = std::make_pair(tmf[54].second, code);
    TI_old = std::visit(my_string(), ftmm.at(tmf[0].first).second);
}
*/

// ----------------------------------------------------------------------------
// Function single_multi_discarded()
// ----------------------------------------------------------------------------
void single_multi_discarded(t_map_str_any & ftmm, t_map_str_any & ftmm_old, int & TI_old, int & AM_after, bool & discard_all)
{
    int code = 3; // single-mapped: 1; multi-mapped: 2; discarded: 3;

    if(TI_old == -1 || TI_old != (std::stoi(std::visit(my_string(), ftmm.at(tmf[0].first).second)))) { // new tag
        discard_all = false;
        int AM_actual = stoi(std::visit(my_string(), ftmm.at(tmf[9].first).second));
        if(AM_actual > AM_after) {
            code = 1;
        }
        else if(AM_actual == AM_after) {
            code = 2;
        }
    }
    else { // same tag
        if(discard_all) {
            code = 3;
        }
        else {
            int AM_actual = std::stoi(std::visit(my_string(), ftmm.at(tmf[9].first).second));
            int AM_old = std::stoi(std::visit(my_string(), ftmm_old.at(tmf[9].first).second));
            if(AM_actual == AM_old) {
                code = std::stoi(std::visit(my_string(), ftmm_old.at(tmf[54].first).second));
            }
            else if(AM_actual < AM_old) {
                code = 3;
                discard_all = true;
            }
        }
    }
    ftmm[tmf[54].first] = std::make_pair(tmf[54].second, code);
    TI_old = std::stoi(std::visit(my_string(), ftmm.at(tmf[0].first).second));
}

// ----------------------------------------------------------------------------
// Function fill_all_tag_mir_fields()
// ----------------------------------------------------------------------------

void fill_and_print_all_tag_mir_fields(t_tag_mir_matched_vv const & mir_tag_out, t_map_str_bool const & org_ids_m,
        options const & options, std::vector<std::ofstream> & tm_files_tab, std::vector<std::ofstream> & tm_files_gff, std::ifstream & multisample_table)
{
    //TODO check that the calculation of PSD field is correct (MAYBE is erroneously overwritten, if so, perhaps the data
    // structure should be reviewed). Correct PSD would be useful (with MSD) to have a more discriminating/selective SMD
    std::string line;
    int line_n = -1;
    t_map_str_any ftmm_old; // Fields Tag Mir (match) Map Old
    bool first_mir = true;
    int AM_after = -1;
    int TI_old = -1;
    bool discard_all = false;
    t_map_str_bool org_bool;

    for(unsigned i = 0; i < mir_tag_out.size(); ++i)
    {
        for(unsigned j = 0; j < mir_tag_out[i].mir_matched.size(); ++j)
        {
            if (mir_tag_out[i].mir_matched[j].mir_prxmir.size() > 0) // this branch is used only if tag is aligned on a prxmir
            {
                for (t_org_prxmir_matched_p const & org_prx: mir_tag_out[i].mir_matched[j].mir_prxmir)
                {
                    if(org_ids_m.count(ALL_ORGANISM) > 0 || org_ids_m.count(org_prx.first) > 0) {
                        for (unsigned w = 0; w < org_prx.second.size(); ++w) {
                            t_map_str_any ftmm; // Fields Tag Mir (match) Map
                            t_prx_mirna_cell const &prxmir = mir_tag_out[i].mir_matched[j].mir->org_prx_m.at(
                                    org_prx.second[w].org)[org_prx.second[w].prxmir_id];
                            t_mir_info_cell mir_info;
                            if (mir_tag_out[i].mir_matched[j].mir->index == prxmir.index5p) {
                                mir_info = prxmir.mirna5p;
                            } else if (mir_tag_out[i].mir_matched[j].mir->index == prxmir.index3p) {
                                mir_info = prxmir.mirna3p;
                            }
                            // Print coordinates in data structure // TODO Decide if remove or keep this line for the tabular file
                            ftmm[tmf[49].first] = std::make_pair(tmf[49].second, i);
                            ftmm[tmf[50].first] = std::make_pair(tmf[50].second, j);
                            ftmm[tmf[6].first] = std::make_pair(tmf[6].second, org_prx.first);
                            ftmm[tmf[51].first] = std::make_pair(tmf[51].second, mir_info.index);
                            ftmm[tmf[52].first] = std::make_pair(tmf[52].second, org_prx.second[w].prxmir_id);
                            // Print info in TagMirMatch data structure
                            fill_tag_mir_match_fields(ftmm, mir_tag_out[i], mir_tag_out[i].mir_matched[j]);
                            //Print miRNA number
                            ftmm[tmf[27].first] = std::make_pair(tmf[27].second,
                                                                 (unsigned) mir_tag_out[i].mir_matched.size());
                            ftmm[tmf[44].first] = std::make_pair(tmf[44].second, (unsigned) org_prx.second.size());
                            // Print difference in alignment score between the j and j+1 (of same organism) aligned mirna
                            if ((j < mir_tag_out[i].mir_matched.size() - 1)) {
                                for(unsigned k = j+1; k < mir_tag_out[i].mir_matched.size(); k++) {
                                    //seqan3::rna15_vector tmp_seq{"TGAGGTAGTAGTTTGTATTATTA"_rna15}; // TI = 19283424
                                    //seqan3::rna15_vector tmp_seq{"AAACAGCACGTAAATATTGGCGA"_rna15}; // TI = 43736
                                    if(mir_tag_out[i].mir_matched[k].mir->org_prx_m.count(org_prx.first) > 0) {
                                        ftmm[tmf[28].first] = std::make_pair(tmf[28].second,
                                                                             mir_tag_out[i].mir_matched[j].align_score -
                                                                             mir_tag_out[i].mir_matched[k].align_score);
                                        break; // first occurrence only
                                    }
                                }

                            }
                            // Print info of aligned mir
                            fill_tag_mir_info_match_fields(ftmm, mir_info,
                                                           mir_tag_out[i].mir_matched[j].mir_tag_size_diff,
                                                           mir_tag_out[i].tag->seq.size());

                            //Print info in TagPremirMatch datastructure
                            ftmm[tmf[45].first] = std::make_pair(tmf[45].second, org_prx.second[w].align_score);
                            ftmm[tmf[46].first] = std::make_pair(tmf[46].second,
                                                                 bool2char(org_prx.second[w].iso5p_canonic));
                            ftmm[tmf[47].first] = std::make_pair(tmf[47].second,
                                                                 bool2char(org_prx.second[w].iso3p_canonic));

                            // Print difference in alignment score between the w and w+1 aligned prxmir
                            if (w < org_prx.second.size() - 1) {
                                ftmm[tmf[48].first] = std::make_pair(tmf[48].second, org_prx.second[w].align_score -
                                                                                     org_prx.second[w + 1].align_score);
                            }
                            // Print prxmir data structure
                            fill_tag_prxmir_info_fields(ftmm, prxmir);
                            //ftmmv.push_back(ftmm);

                            // Actual print to file
                            //if(org_prx.second.size() > 1 && w<(org_prx.second.size()-1)) {
                                //if(org_prx.second[w].align_score == org_prx.second[w + 1].align_score) {
                                    //same_mir_diff_prx = true;
                                //}
                            //}
                            if(w+1 < org_prx.second.size()) {
                                AM_after = std::stoi(std::visit(my_string(), ftmm.at(tmf[9].first).second));
                            }
                            else if(j+1 < mir_tag_out[i].mir_matched.size()) {
                                AM_after = mir_tag_out[i].mir_matched[j+1].align_score;
                            }
                            else {
                                AM_after = -1;
                            }
                            single_multi_discarded(ftmm, ftmm_old, TI_old, AM_after, discard_all);
                            print_tag_mir_tab(tm_files_tab, ftmm, multisample_table, line, line_n, options);

                            if( j > 0 && ftmm.at(tmf[0].first).second == ftmm_old.at(tmf[0].first).second) { // still same tag
                                first_mir = false;
                                if(ftmm.at(tmf[6].first).second != ftmm_old.at(tmf[6].first).second) { // if different organism
                                    if(org_bool.count(std::visit(my_string(), ftmm.at(tmf[6].first).second)) <= 0) { // if organism not yet considered
                                        first_mir = true;
                                        org_bool[std::visit(my_string(), ftmm.at(tmf[6].first).second)] = true; // record organism checked (useful for FILTER=)
                                    }
                                }
                                else if(ftmm.at(tmf[30].first).second == ftmm_old.at(tmf[30].first).second) { // if organism already considered but same miRNA (same miRNA but different prxmir)
                                    first_mir = true;
                                }
                            }
                            else { // new tag considered
                                first_mir = true;
                                org_bool.clear();
                                org_bool[std::visit(my_string(), ftmm.at(tmf[6].first).second)] = true; // record organism checked (useful for FILTER=)
                            }
                            print_tag_mir_gff(tm_files_gff, ftmm, multisample_table, line, first_mir, options);
                            ftmm_old = ftmm;
                        }
                    }
                }
            } else // this branch is used if tag is not aligned on a prxmir because it is shorter
            {
                for (t_org_prx_p const & org_prx: mir_tag_out[i].mir_matched[j].mir->org_prx_m)
                {
                    if(org_prx.second.size() > 0) // if prxmir exist (maybe this check is useless)
                    {
                        if(org_ids_m.count(ALL_ORGANISM) > 0 || org_ids_m.count(org_prx.first) > 0) {
                            for (unsigned l = 0; l < org_prx.second.size(); ++l) {
                                t_map_str_any ftmm;
                                t_prx_mirna_cell const &prxmir = org_prx.second[l];
                                t_mir_info_cell mir_info;
                                if (mir_tag_out[i].mir_matched[j].mir->index == prxmir.index5p) {
                                    mir_info = prxmir.mirna5p;
                                } else if (mir_tag_out[i].mir_matched[j].mir->index == prxmir.index3p) {
                                    mir_info = prxmir.mirna3p;
                                }
                                // Print coordinate in data structure //TODO Decide if remove or keep this line for the tabular file
                                ftmm[tmf[49].first] = std::make_pair(tmf[49].second, i);
                                ftmm[tmf[50].first] = std::make_pair(tmf[50].second, j);
                                ftmm[tmf[6].first] = std::make_pair(tmf[6].second, org_prx.first);
                                ftmm[tmf[51].first] = std::make_pair(tmf[51].second, mir_info.index);
                                ftmm[tmf[52].first] = std::make_pair(tmf[52].second, l);
                                // Print info in TagMirMatch data structure
                                fill_tag_mir_match_fields(ftmm, mir_tag_out[i], mir_tag_out[i].mir_matched[j]);
                                // Print miRNA number
                                ftmm[tmf[27].first] = std::make_pair(tmf[27].second,
                                                                     (unsigned) mir_tag_out[i].mir_matched.size());
                                // Print difference in alignment score between the j and j+1 (of same organism) aligned mirna
                                if ((j < mir_tag_out[i].mir_matched.size() - 1)) {
                                    for(unsigned k = j+1; k < mir_tag_out[i].mir_matched.size(); k++) {
                                        if(mir_tag_out[i].mir_matched[k].mir->org_prx_m.count(org_prx.first) > 0) {
                                            ftmm[tmf[28].first] = std::make_pair(tmf[28].second,
                                                                                 mir_tag_out[i].mir_matched[j].align_score -
                                                                                 mir_tag_out[i].mir_matched[k].align_score);
                                            break; // first occurrence only
                                        }
                                    }

                                }
                                // Print info of aligned mir
                                fill_tag_mir_info_match_fields(ftmm, mir_info,
                                                               mir_tag_out[i].mir_matched[j].mir_tag_size_diff,
                                                               mir_tag_out[i].tag->seq.size());
                                //Print premir datastructure
                                fill_tag_prxmir_info_fields(ftmm, prxmir);
                                //ftmmv.push_back(ftmm);

                                // Actual print to file
                                if(l+1 < org_prx.second.size()) {
                                    AM_after = std::stoi(std::visit(my_string(), ftmm.at(tmf[9].first).second));
                                }
                                else if(j+1 < mir_tag_out[i].mir_matched.size()) {
                                    AM_after = mir_tag_out[i].mir_matched[j+1].align_score;
                                }
                                else {
                                    AM_after = -1;
                                }
                                single_multi_discarded(ftmm, ftmm_old, TI_old, AM_after, discard_all);
                                print_tag_mir_tab(tm_files_tab, ftmm, multisample_table, line, line_n, options);

                                if( j > 0 && ftmm.at(tmf[0].first).second == ftmm_old.at(tmf[0].first).second) { // still same tag
                                    first_mir = false;
                                    if(ftmm.at(tmf[6].first).second != ftmm_old.at(tmf[6].first).second) { // if different organism
                                        if(org_bool.count(std::visit(my_string(), ftmm.at(tmf[6].first).second)) <= 0) { // if organism not yet considered
                                            first_mir = true;
                                            org_bool[std::visit(my_string(), ftmm.at(tmf[6].first).second)] = true; // record organism checked (useful for FILTER=)
                                        }
                                    }
                                    else if(ftmm.at(tmf[30].first).second == ftmm_old.at(tmf[30].first).second) { // if organism already considered but same miRNA (same miRNA but different prxmir)
                                        first_mir = true;
                                    }
                                }
                                else { // new tag considered
                                    first_mir = true;
                                    org_bool.clear();
                                    org_bool[std::visit(my_string(), ftmm.at(tmf[6].first).second)] = true; // record organism checked (useful for FILTER=)
                                }
                                print_tag_mir_gff(tm_files_gff, ftmm, multisample_table, line, first_mir, options);
                                ftmm_old = ftmm;
                            }
                        }
                    }
                }
            }
        }
    }
}

void print_discarded_tags(t_tag & tag, options const & options)
{
    _V(options, "STEP: Start printing IsomiR-SEA discarded tags");
    int index_id = 0;
    t_tag_v tags_out;

    for(unsigned i = 0; i < tag.tag_v.size(); ++i) {
        if(tag.tag_v[i].mir_matched_id == -1)
        {
            tag.tag_v[i].reason_discard = "NOTALL";
            tags_out.push_back(tag.tag_v[i]);
        }
    }
    for(unsigned i = 0; i < tag.discarded_tag_v.size(); ++i) {
        tag.discarded_tag_v[i].reason_discard = "LENGTH";
        tags_out.push_back(tag.discarded_tag_v[i]);
    }
    std::sort(tags_out.begin(), tags_out.end(), [](t_tag_cell a, t_tag_cell b) { return a.seq < b.seq; });
    sequence_file_output fout{options.path_discarded_tags};
    for(unsigned i=0; i<tags_out.size(); ++i) {
        std::vector<dna15> seq;
        for(unsigned k=0; k<tags_out[i].seq.size(); ++k) { // convert from rna15 to dna15 for YARA computation
            dna15 nucl;
            nucl.assign_rank(tags_out[i].seq[k].to_rank()); // rank is alphabet_rank_t<seqan3::dna15>
            seq.push_back(nucl);
        }
        //std::string header = "ID:" + std::to_string(index_id) + "|" +
        std::string header = "ID:" + std::to_string(tags_out[i].index) + "|" +
                             "CN:" + std::to_string(tags_out[i].count) + "|" +
                             "RS:" + tags_out[i].reason_discard;
        fout.emplace_back(seq, header, tags_out[i].qual);
        ++index_id;
    }
}

#endif //ISOMIR_SEA_OUTPUT_H
