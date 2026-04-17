/*
 * =============================================================================
 * main.cpp
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef SEQAN3_ISOMIR_SEA
#define SEQAN3_ISOMIR_SEA

// ============================================================================
// Macro utility
// ============================================================================

#define _V(_opt, _str) {if(_opt.verbose>0) std::cerr << _str << "\n";}
#define _VV(_opt, _str) {if(_opt.verbose>1) std::cerr << _str << "\n";}
#define _VVV(_opt, _str) {if(_opt.verbose>2) std::cerr << _str << "\n";}

// ----------------------------------------------------------------------------
// Utility headers
// ----------------------------------------------------------------------------

#include <boost/algorithm/string.hpp>
//#include <boost/variant/static_visitor.hpp>
#include <range/v3/range/conversion.hpp>
#include <omp.h>
#include <regex>

// ----------------------------------------------------------------------------
// SeqAn3 headers
// ----------------------------------------------------------------------------

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

// ----------------------------------------------------------------------------
// SeqAn2 headers
// ----------------------------------------------------------------------------

#include <seqan/find.h>

// ----------------------------------------------------------------------------
// Experimental headers
// ----------------------------------------------------------------------------

#include <experimental/filesystem> // C++14

// ----------------------------------------------------------------------------
// Cereal headers
// ----------------------------------------------------------------------------

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/variant.hpp>
#include <cereal/types/utility.hpp>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "mirna.h"
#include "prx_mirna.h"
#include "seed_mirna.h"
#include "tag.h"
#include "prxmir_matched.h"
#include "mir_matched.h"
#include "tag_mirna.h"
#include "options.h"
#include "input.h"
#include "iupac_map.h"
#include "core.h"
#include "output.h"
#include "stopwatch.h"

using namespace seqan3;
namespace fs = std::experimental::filesystem;
typedef std::experimental::filesystem::path path;

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main(int argc, char **argv)
{
    // Argument Parser
    if (argc == 1)
    {
        std::cerr << "Type " << argv[0] << " -h to get the parameters table" << "\n";
        return 1;
    }
    argument_parser parser("isomiR-SEA", argc, argv);
    options options{};
    setup_argument_parser(parser, options); // prepare the option structure
    try
    {
        parser.parse(); // trigger parsing
    }
    catch (parser_invalid_argument const & ext)
    {
        std::cerr << "Exception caught: " << ext.what() << "\n";
        return 1;
    }
    if(!check_parse_result(options)) {
        return 1;
    }
    _V(options, "STEP: Verify paths and create out dir");
    if(check_out_path(options)) {
        return -1;
    }
    stopwatch<> total_time("Total Time");

    _V(options, "INFO: Input file selected by user is: " << options.in_file_tags);
    _V(options, "INFO: Label of out directory is: " << options.out_file_label);
    _V(options, "INFO: Out directory selected by user is: " << options.path_out_files);

    path dir = path(options.in_file_tags).parent_path();
    path out_dir(options.path_out_files);
    out_dir /= path(options.out_file_label);
    fs::create_directory(out_dir);

    // Check if paths are absolute
    if(!dir.is_absolute()) {
        dir = fs::current_path() / dir;
    }
    if(!out_dir.is_absolute()) {
        out_dir = fs::current_path() / out_dir;
    }
    _V(options, "INFO: Input directory is: " << dir.generic_string());
    _V(options, "INFO: Out directory is: " << out_dir.generic_string());

    fs::path p_log = out_dir / "align.log";
    std::ofstream log_file;
    log_file.open(p_log);
    print_options(log_file, options);  // print options used for isomiR-SEA run on log file
    log_file << "\n#Alignment features detected during the alignment of a Tag over miRs" << std::endl;
    for(unsigned i = 0; i < tmf.size(); ++i) // print alignment features detected during the alignment on log file
    {
        log_file << i << "\t" << tmf[i].first << "\t" << tmf[i].second << std::endl;
    }

    _V(options, "STEP: Read reference files");
    stopwatch<> in_load_time("Load time of input files");
    t_seed_v seeds;
    t_org_prxmir_m ref_prxmir_db;
    t_map_str_bool org_ids_m;
    extract_elements(org_ids_m, options.specie_codes, options);

    if(options.path_load_serialized.empty())
    {
        load_reference_db(seeds, ref_prxmir_db, org_ids_m, options);
        if(!options.path_store_serialized.empty())
        {
            _V(options, "STEP: Serializing and storing reference files");
            fs::path p_db_cereal = options.path_store_serialized;
            std::ofstream db_cereal_file_out;
            db_cereal_file_out.open(p_db_cereal, std::ios_base::binary);
            cereal::BinaryOutputArchive oarchive(db_cereal_file_out);
            oarchive(ref_prxmir_db, seeds, MAX_MIR_SIZE);
        }
    }
    else
    {
        _V(options, "STEP: Loading serialized reference files");
        fs::path p_db_cereal = options.path_load_serialized;
        std::ifstream db_cereal_file_in;
        db_cereal_file_in.open(p_db_cereal, std::ios_base::binary);
        cereal::BinaryInputArchive iarchive(db_cereal_file_in);
        iarchive(ref_prxmir_db, seeds, MAX_MIR_SIZE);
    }

    // Remove miRNA-sequence and seeds according to user-selected organisms (in order to save lot of alignment time)
    prune_datastructure(seeds, org_ids_m);

    // Verify that the inserted specie_codes are available in the dataset
    if(org_ids_m.count(ALL_ORGANISM) <= 0) {
        for (t_pair_str_bool const &org_bool : org_ids_m) {
            if (ref_prxmir_db[org_bool.first].empty()) {
                std::cerr << "User-requested organism '" << org_bool.first
                             << "' is not present in the referenceDB, exiting...";
                return 1;
            }
        }
    }
    _V(options, "INFO: unique seeds " << seeds.size());
    if(options.verbose > 3) // print a file storing all the reference for selected organisms
    {
        fs::path p_mir_db = out_dir / "mir_db.txt";
        std::ofstream mir_db_file;
        mir_db_file.open(p_mir_db);
        print_seed_mir_prx(mir_db_file, seeds);
    }

    _V(options, "STEP: Read Tag file");
    t_tag tag;
    if(fs::path(options.in_file_tags).extension() == ".tag" || fs::path(options.in_file_tags).extension() == ".tagq") {
        load_tags(tag, options);
    }
    else {
        // FASTA or FASTQ — check first line to detect format
        std::ifstream infile(options.in_file_tags);
        std::string line;
        // bool is_bioseqzip = false;
        bool is_mirtrace = false;
        std::regex pattern("_x\\d+ ");

        if (infile && std::getline(infile, line)) {
            if (!line.empty() && (line[0] == '>' || line[0] == '@')) {
                // Check if the ID line contains '|'
                //is_bioseqzip = line.find('|') != std::string::npos;
                is_mirtrace = std::regex_search(line, pattern);
            }
        }
        infile.close();

        if (is_mirtrace) {
            if (options.mirtrace_input == false) {
                std::cout << "Warning: Input file seems to be in miRTrace format, but --mirtrace-input option was not set. Proceeding with miRTrace format." << std::endl;
            }
            load_mirtrace_fastx(tag, options);
        } else {
            load_tags_fastx(tag, options);
        }

        /*
        if (!is_bioseqzip && options.mirtrace_input) {
            load_mirtrace_fastx(tag, options);
            
        } else {
            load_tags_fastx(tag, options);
        */

        /*
            else{
                debug_stream << "Error: Input file format not recognized. Please provide a .tag, .tagq, .fasta or .fastq file.\n";
                return 1;
            }
            */
    }
    in_load_time.print(log_file);
    _V(options, "INFO: Tags used for IsomiR-SEA alignment = " << tag.tag_v.size());
    _V(options, "INFO: Tags discarded = " << tag.discarded_tag_v.size());
    log_file << "INFO: Tags used for IsomiR-SEA alignment = " << tag.tag_v.size() << std::endl;
    log_file << "INFO: Tags discarded = " << tag.discarded_tag_v.size() << std::endl;

    _V(options, "STEP: Start computation of IsomiR-SEA alignment");
    stopwatch<> isea_time("IsomiR-SEA alignment Time");
    int count_align_mir = 0;
    t_tag_mir_matched_vv mir_tag_out;
    std::vector<t_tag_mir_matched_vv> mir_buffers;
    std::vector<int> count_align_mir_partial;

    omp_set_num_threads(options.threads);

    #pragma omp parallel
    {
        auto nthreads = options.threads;
        auto id = omp_get_thread_num();

        // Correctly set the number of buffers
        #pragma omp single
        {
            mir_buffers.resize(nthreads);
            count_align_mir_partial.resize(nthreads, 0);
        };
        // Perform scan
        #pragma omp for schedule(static)
        for(unsigned i = 0; i < tag.tag_v.size(); ++i)
        {
            t_tag_mir_matched_v tag_mir_matched_v;
            //debug_stream << i << "\t" << tag.tag_v[i].seq << "\n";
            scan_mir(tag_mir_matched_v.mir_matched, tag.tag_v[i], seeds, org_ids_m, options);
            if(tag_mir_matched_v.mir_matched.size() > 0)
            {
                tag_mir_matched_v.tag = & tag.tag_v[i];
                mir_buffers[id].push_back(tag_mir_matched_v);
                tag.tag_v[i].mir_matched_id = count_align_mir_partial[id]; //TODO This id is not unique in multithread computation (a possible solution could be a new for over the mir_tag_out size)
                count_align_mir_partial[id] += tag_mir_matched_v.mir_matched.size();
            }
        }
        // Combine buffers together
        #pragma omp single
        {
            for(auto & buffer : mir_buffers)
            {
                move(buffer.begin(), buffer.end(), back_inserter(mir_tag_out));
            }
            for(unsigned i = 0; i < count_align_mir_partial.size(); ++i)
            {
                count_align_mir += count_align_mir_partial[i];
            }
        }
    };
    for(unsigned  i = 0; i < mir_tag_out.size(); ++i)
    {
        mir_tag_out[i].tag->mir_matched_id = i;
    }
    isea_time.print(log_file); // print time required to performs all the alignments of selected tags over all miRNAs
    _V(options, "STEP: End computation of IsomiR-SEA alignment");
    _V(options, "INFO: Num Total Tag Seqs = " << tag.tag_v.size() + tag.discarded_tag_v.size());
    _V(options, "INFO: Num Total Reads Seqs = " << tag.read_count + tag.invalid_read_count);
    _V(options, "INFO: Num Align Tag Seqs = " << mir_tag_out.size());
    _V(options, "INFO: Num Align Mir Seqs = " << count_align_mir);
    log_file << "INFO: Num Total Tag Seqs = " << tag.tag_v.size() + tag.discarded_tag_v.size() << "\n";
    log_file << "INFO: Num Total Reads Seqs = " << tag.read_count + tag.invalid_read_count << "\n";
    log_file << "INFO: Num Align Tag Seqs = " << mir_tag_out.size() << "\n";
    log_file << "INFO: Num Align Mir Seqs = " << count_align_mir << "\n";

    _V(options, "STEP: Start printing IsomiR-SEA alignments");

    t_map_str_bool file_types_m;
    extract_elements(file_types_m, options.out_formats, options); // this line create a list of file formats to be printed
    //for(auto & file_type: file_types_m) { debug_stream << file_type.first << " = " << file_type.second << "\n"; }
    // TODO eventually add possibility to print other output file formats (user-selection)

    // Preparing output file descriptors
    std::ifstream multisample_table;
    std::vector<fs::path> tm_all_tab_v, tm_all_gff_v;
    std::string line;
    if(options.multi_sample_tab.empty())
    {
        fs::path tm_all_tab = out_dir / "tag_mir-all.tab";
        tm_all_tab_v.push_back(tm_all_tab);
        fs::path tm_all_gff = out_dir / "tag_mir-all.gff";
        tm_all_gff_v.push_back(tm_all_gff);
    }
    else
    { // Retrieve filenames from first line of BioSeqZip output tabular file (read-count per sample)
        t_vect_str filenames;
        multisample_table.open(options.multi_sample_tab);
        std::getline(multisample_table, line);
        boost::split(filenames, line, boost::is_any_of("\t"));
        for(unsigned i=1; i<filenames.size(); ++i) {
            std::string filename_tab = filenames[i] + ".tab";
            std::string filename_gff = filenames[i] + ".gff";
            fs::path path_tab = out_dir / filename_tab;
            tm_all_tab_v.push_back(path_tab);
            fs::path path_gff = out_dir / filename_gff;
            tm_all_gff_v.push_back(path_gff);
        }
    }
    // Open tab file(s) and print header
    std::vector<std::ofstream> tm_files_tab;
    tm_files_tab.resize(tm_all_tab_v.size());
    for(unsigned i=0; i<tm_all_tab_v.size(); ++i) {
        tm_files_tab[i].open(tm_all_tab_v[i]);
        for(unsigned j = 0; j < tmf.size(); ++j) {
            tm_files_tab[i] << tmf[j].second << "\t";
        }
        tm_files_tab[i] << "\n";
    }
    // Open gff file(s)
    std::vector<std::ofstream> tm_files_gff;
    tm_files_gff.resize(tm_all_gff_v.size());
    for(unsigned i=0; i<tm_all_gff_v.size(); ++i) {
        tm_files_gff[i].open(tm_all_gff_v[i]);
    }

    //t_vect_map_str_any fields_tag_mir_v;
    fill_and_print_all_tag_mir_fields(mir_tag_out, org_ids_m, options, tm_files_tab, tm_files_gff, multisample_table);

    // Close all file descriptiors
    for(unsigned i=0; i<tm_all_tab_v.size(); ++i) {
        tm_files_tab[i].close();
    }
    for(unsigned i=0; i<tm_all_gff_v.size(); ++i) {
        tm_files_gff[i].close();
    }
    multisample_table.close();
    //fields_tag_mir_v.clear();

    _V(options, "STEP: End printing IsomiR-SEA alignments");

    // Print discarded tags
    if(!options.path_discarded_tags.empty()) {
        print_discarded_tags(tag, options);
    }

    total_time.print(log_file);

    return 0;
}

#endif // SEQAN3_ISOMIR_SEA
