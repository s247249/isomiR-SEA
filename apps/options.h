/*
 * =============================================================================
 * options.h
 *
 * Options and configuration parameters for the tool.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_OPTIONS_H
#define ISOMIR_SEA_OPTIONS_H

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

//#include <iostream>
//#include <string>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan3/argument_parser/argument_parser.hpp>

// Returns a temporary filename.
#ifndef SEQAN_TEMP_FILENAME
#define SEQAN_TEMP_FILENAME() (::temp_file_name())
#endif

using namespace seqan3;

// ----------------------------------------------------------------------------
// Class options
// ----------------------------------------------------------------------------
struct options
{
// Name of miR database
    std::string mir_db_name{};
// Version of miR database
    std::string mir_db_version{};
// Website of miR database
    std::string mir_db_web{};
// Name of input file storing Tags (full path or searched in the current folder)
    std::string in_file_tags{};
// If set, the input file is in miRTrace format (4 columns: tag_id, sequence, count, and miRNA annotation)
    bool mirtrace_input{false};
// Name of input file storing mature/star miRs (if passed alone both mature/star are acquired and will be tryed to discriminate the type (full path or searched in the current folder)
    std::string in_file_mature{};
// Name of input file storing star miRs that if passed separately can be readed and included in the mature structure (full path or searched in the current folder)
    std::string in_file_star{};
// Name of input file storing GFF3 with miR/premiR annotated on ref genome (full path or searched in the current folder)
    std::string in_file_gff{};
// Type field of miRNA in GFF3 input file
    std::string type_mir_gff{"miRNA"};
// Name of input file storing BED with miR/premiR annotated on ref genome (full path or searched in the current folder)
    //std::string in_file_bed{};
// Name of input file storing premiRs (full path or searched in the current folder)
    //std::string in_file_premir{};
// Load referenceDB from serialized file
    std::string path_load_serialized{};
// Serialize and store referenceDB
    std::string path_store_serialized{};
// Name of input file storing primiRs (full path or searched in the current folder)
    std::string in_file_primir{};
// Unformatted name of mature miRNAs (format is ORG-MIRNAME and eventually Prologuos, 5p/3p, and *
    bool unform_mature_names{false};
// If this parameter is True the conversion between DNA to RNA is avoided
    bool dna_or_rna{true};
// Path to a tabular file containing #reads per sample (which means <in-file-tag> file is multi-sample)
    std::string multi_sample_tab{};
// Name of common part of the output file (Default will be used the last folder of the path storing the dataset)
    std::string out_file_label{};
// Formats of the output file (tab and gff are available)
    std::string out_formats{"tab-gff"};
// Path where to generate the folders storing output files (If empty will be used the database name)
    std::string path_out_files{};
// Path where to save discarded tag to be used in YARA
    std::string path_discarded_tags{};
// Temporary directory where to save intermediate files. Default: use \tmp of OS.
    std::string tmp_dir{};

//Allign unmapped tag on premir
    //bool dp_allign_on_premirs;
// Gap scheme to be used in the N-W alignment (default: affine(0))
    //unsigned affine_linear_dgs{0u};
// Gap open and extend costs
    //double gap_open{-6.0};
    //double gap_extend{-2.0};
// Score for a match
    //double match_score{4.0};
// Score for a mismatch
    //double mismatch_score{-2.0};
// Use the global local or global-Unconstrained algorithm (default: global(0) - local(1) )
    //bool global_local{false};
// type used for the global-unconstrained alignment AlignConfig <TTop, Tleft, TRight, TDown>
    //bool un_top{false};
    //bool un_left{false};
    //bool un_right{false};
    //bool un_down{false};

// Specie codes to be evaluated during the alignment at the same time
    std::string specie_codes{}; // only alphanums can be used to represents the organism code, the other chars are considered separators
// Minimum size of ungapped alignment, starting from the seed, extending the alignment
    unsigned min_size_aln_stp1{10u};
// Minimum alignment score for considering a tag expression of a miR
    int min_align_score{7};
// Start position of the seed
    int seed_start{1u};
// End position of the seed
    int seed_end{6u};
// Max index in tag position for starting the seed alignment
    unsigned max_start_pos_tag{5u};
// Minimum size of tag to be considered for the alignment
    unsigned min_size_tag{15u};
// Threshold used to select select high quality multimapped tags
    int thr_select_tags{12};
// Number of mismatches allowed between miRNA and tags
    int mismatches_out_seed{3};
// Number of mismatches allowed between miRNA seed and tags
    int mismatches_in_seed{0};

// Time used for an hard timeout
    int time_limit{-1};
// verbose(0) no outputs,
// verbose(1) Displays global statistics,
// verbose(2) Displays extensive statistics for each batch of reads,
// verbose(3) Debug output.
    unsigned verbose{0u};
#ifdef _OPENMP
    // number of threads forced
    unsigned threads{std::thread::hardware_concurrency()}; // omp_get_num_threads()
    // number of threads detected
    unsigned threads_count{std::thread::hardware_concurrency()};
#else
    // number of threads forced
    unsigned threads{1u};
    // number of threads detected
    unsigned threads_count{1u};
#endif
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function setup_argument_parser()
// ----------------------------------------------------------------------------
void setup_argument_parser(argument_parser & parser, options & options)
{
    parser.info.app_name = "isomiR-SEA";
    parser.info.short_description = "an RNA-Seq analysis tool for miRNAs/isomiRs expression level profiling and miRNA-mRNA interaction sites evaluation";
    // TODO: SeqAn3 counterpart for SeqAn2 setCategory?
    parser.info.version = "2.0";
    parser.info.date = "2019";
    parser.info.synopsis = {"./isomir_sea <\\fI--in-file-tags inFileTags --in-file-primir inFilePrimir --in-file-gff inFileGff\\fP> [\\fI--out-file-label outFileLabel\\fP] [\\fI-parameters\\fP]"}; // old project was: "./isomir_sea <\\fI-ift inFileTags -ifm inFileMature\\fP> [\\fI-ofl outFileLabel\\fP] [\\fI -parameters\\fP]";
    parser.add_option(options.verbose, 'V' /*"v"*/, "verbose", "verbose(0) no outputs, verbose(1) Displays global statistics, verbose(2) Displays extensive statistics for each batch of reads, verbose(3) Debug output.");
    parser.add_section("isomiR-SEA Alignment Options");
    // Second parameter (short_name) can only be a single char in seqan3, so I simply used incremental characters
    //parser.add_option(options.dp_allign_on_premirs, 'a' /*dap*/, "dp-allign-on-premirs", "Allign unmapped tag on premir.");
    //parser.add_option(options.affine_linear_dgs, 'b' /*ald*/, "affine-linear-dgs", "Chose the gap scheme affine(0) linear(1) or dynamic(2) to be used in the alignment.");
    //parser.add_option(options.gap_open, 'c' /*"go"*/, "gap-open", "Gap open cost for the Smith-Waterman alignment.");
    //parser.add_option(options.gap_extend, 'd' /*"ge"*/, "gap-extend", "Gap extend cost for the Smith-Waterman alignment.");
    //parser.add_option(options.match_score, 'e' /*"ma"*/, "match-score", "Score of a match for the Smith-Waterman alignment.");
    //parser.add_option(options.mismatch_score, 'f' /*"mm"*/, "mismatch-score", "Score of a mismatch for the Smith-Waterman alignment.");
    //parser.add_option(options.global_local, 'g' /*"gl"*/, "global-local", "Use the global local(0) or global-unconstrained(1) algorithm.");
    //parser.add_option(options.un_top, 'H' /*"ut"*/, "un-top", "type used for the global-unconstrained alignment AlignConfig TTop."); // -h is already used by argument_parser itself
    //parser.add_option(options.un_left, 'i' /*"ul"*/, "un-left", "type used for the global-unconstrained alignment AlignConfig TLeft.");
    //parser.add_option(options.un_right, 'l' /*"ur"*/, "un-right", "type used for the global-unconstrained alignment AlignConfig TRight.");
    //parser.add_option(options.un_down, 'm' /*"ud"*/, "un-down", "type used for the global-unconstrained alignment AlignConfig TDown.");
    parser.add_option(options.specie_codes, 'n' /*"sc"*/, "specie-codes", "Specie codes to be evaluated during the alignment at the same time.", option_spec::required);
    parser.add_option(options.min_size_aln_stp1, 'o' /*"msas"*/, "min-size-aln-stp1", "Minimum size of ungapped alignment, starting from the seed, extending the alignment.");
    parser.add_option(options.min_align_score, 'p' /*"mas"*/, "min-align-score", "Minimum alignment score for considering a tag expression of a miR.");
    parser.add_option(options.seed_start, 'q' /*"ss"*/, "seed-start", "Start position of the seed.");
    parser.add_option(options.seed_end, 'r' /*"se"*/, "seed-end", "End position of the seed.");
    parser.add_option(options.max_start_pos_tag, 's' /*"mspt"*/, "max-start-pos-tag", "Max index in tag position for starting the seed alignment.");
    parser.add_option(options.min_size_tag, 't' /*"mst"*/, "min-size-tag", "Minimum size of tag to be considered for the alignment.");
    parser.add_option(options.thr_select_tags, 'u' /*"tst"*/, "thr-select-tags", "Threshold used to select select high quality multimapped tags.");
    parser.add_option(options.mismatches_out_seed, 'v' /*"mos"*/, "mismatches-out-seed", "Number of mismatches allowed between miRNA and tags.");
    parser.add_option(options.mismatches_in_seed, 'z' /*"mis"*/, "mismatches-in-seed", "Number of mismatches allowed between miRNA seed and tags.");

    parser.add_section("Input Options");
    parser.add_option(options.in_file_tags, '1' /*"ift"*/, "in-file-tags", "Name of input file storing Tags (full path or searched in the current folder).", option_spec::required);
    parser.add_flag(options.mirtrace_input, 'm' /*"mt"*/, "mirtrace-input", "If set, the input file is in miRTrace format (4 columns: tag_id, sequence, count, and miRNA annotation).");
    parser.add_option(options.in_file_mature, '2' /*"ifm"*/, "in-file-mature", "Name of input file storing mature/star miRs (if passed alone both mature/star are acquired and will be tryed to discriminate the type (full path or searched in the current folder).");
    parser.add_option(options.in_file_star, '3' /*"ifs"*/, "in-file-star", "Name of input file storing star miRs that if passed separately can be readed and included in the mature structure (full path or searched in the current folder) [OPTIONAL].");
    parser.add_option(options.in_file_gff, '4' /*"ifg"*/, "in-file-gff", "Name of input file storing GFF3 with miR/premiR annotated on ref genome (full path or searched in the current folder) [OPTIONAL].");
    parser.add_option(options.type_mir_gff, '5' /*"tmg"*/, "type-mir-gff", "Type field of miRNA in GFF3 input file [OPTIONAL].");
    //parser.add_option(options.in_file_bed, '6' /*"ifb"*/, "in-file-bed", "Name of input file storing BED with miR/premiR annotated on ref genome (full path or searched in the current folder) [OPTIONAL].");
    //parser.add_option(options.in_file_premir, '7' /*"ifpe"*/, "in-file-premir", "Name of input file storing Pre-miRs (full path or searched in the current folder) [OPTIONAL].");
    parser.add_option(options.path_load_serialized, '6', "load-serialized", "Load referenceDB from serialized file [OPTIONAL].");
    parser.add_option(options.path_store_serialized, '7', "store-serialized", "Serialize and store referenceDB [OPTIONAL].");
    parser.add_option(options.in_file_primir, '8' /*"ifpi"*/, "in-file-primir", "Name of input file storing Pri-miRs (full path or searched in the current folder) [OPTIONAL].");
    parser.add_option(options.dna_or_rna, '9' /*"dor"*/, "dna-or-rna", "If this parameter is True the conversion between DNA to RNA is avoided.");
    parser.add_option(options.unform_mature_names, '0' /*"umn"*/, "unform-mature-names", "Unformatted name of mature miRNAs (format is ORG-MIRNAME and eventually Prologuos, 5p/3p, and *.");
    parser.add_option(options.multi_sample_tab, '@', "multi-sample-tab", "Path to a tabular file containing #reads per sample (which means <in-file-tag> file is multi-sample)");

    parser.add_section("Output Options");
    parser.add_option(options.out_file_label, 'A' /*"ofl"*/, "out-file-label", "Name of common part of the output file (Default will be used the last folder of the path storing the dataset).");
    parser.add_option(options.out_formats, 'B' /*"of"*/, "out-formats", "Formats of the output file (tab and gff are available).");
    parser.add_option(options.path_out_files, 'C' /*"pof"*/, "path-out-files", "Path were to generate the folders storing output files (If empty will be used the current folder).");
    parser.add_option(options.tmp_dir, 'D' /*"td"*/, "tmp-dir", "Temporary directory where to save intermediate files. Default: use \\tmp of OS.");
    parser.add_option(options.path_discarded_tags, 'E', "path-discarded-tags", "Path where to save discarded tag to be used in YARA");

    // Setup miRNA DB info.
    parser.add_section("miRNA database Info");
    parser.add_option(options.mir_db_name, 'F' /*"dbn"*/, "mir-db-name", "Name of miR database.");
    parser.add_option(options.mir_db_version, 'G' /*"dbv"*/, "mir-db-version", "Version of miR database.");
    parser.add_option(options.mir_db_web, 'I' /*"dbw"*/, "mir-db-web", "Website of miR database.");

    // Setup performance options.
    parser.add_section("Performance Options");
    parser.add_option(options.time_limit, 'L' /*"tl"*/, "time-limit", "Time limit for the program.");
#ifdef _OPENMP
    parser.add_option(options.threads, 'M' /*"t"*/, "threads", "Specify the number of threads to use.",  option_spec::standard, arithmetic_range_validator{1, (double)std::thread::hardware_concurrency() + 1});
#else
    parser.add_option(options.threads, 'M' /*"t"*/, "threads", "Specify the number of threads to use.",  option_spec::standard, arithmetic_range_validator{1, 1});
#endif
}

// ----------------------------------------------------------------------------
// Function temp_file_names() imported from SeqAn2, see: /include/seqan/basic/debug_test_system.h
// ----------------------------------------------------------------------------

static::std::vector<std::string> & temp_file_names()
{
    static::std::vector<std::string> filenames;
    return filenames;
}

// ----------------------------------------------------------------------------
// Function temp_file_name() imported from SeqAn2, see: /include/seqan/basic/debug_test_system.h
// Return the path to a temporary file, in a static buffer in this function. This is not thread safe!
// ----------------------------------------------------------------------------

inline
const char * temp_file_name()
{
    static char file_name_buffer[1000];
#ifdef STDLIB_VS
    static char file_path_buffer[1000];
    //  Gets the temp path env string (no guarantee it's a valid path).
    DWORD dw_ret_val = 0;
    dw_ret_val = GetTempPath(1000,            // length of the buffer
                           file_path_buffer); // buffer for path
    if (dw_ret_val > 1000 || (dw_ret_val == 0))
    {
        std::cerr << "GetTempPath failed" << std::endl;
        exit(1);
    }

    UINT u_ret_val   = 0;
    u_ret_val = GetTempFileName(file_path_buffer,   // directory for tmp files
                              TEXT("SEQAN."),   // temp file name prefix
                              0,                // create unique name
                              file_name_buffer);  // buffer for name

    if (u_ret_val == 0)
    {
        std::cerr << "GetTempFileName failed" << std::endl;
        exit(1);
    }

    DeleteFile(file_name_buffer);
    CreateDirectoryA(file_name_buffer, NULL);
    temp_file_names().push_back(file_name_buffer);
    strcat(file_name_buffer, "\\test_file");
    return file_name_buffer;

#else  // ifdef STDLIB_VS
    strcpy(file_name_buffer, "/tmp/SEQAN.XXXXXXXXXXXXXXXXXXXX");
    mode_t cur_umask = umask(S_IRWXO | S_IRWXG);  // to silence Coverity warning
    int _tmp = mkstemp(file_name_buffer);
    (void) _tmp;
    umask(cur_umask);
    unlink(file_name_buffer);
    mkdir(file_name_buffer, 0777);

    temp_file_names().push_back(file_name_buffer);

    strcat(file_name_buffer, "/test_file");
    return file_name_buffer;

#endif  // ifdef STDLIB_VS
}


// ----------------------------------------------------------------------------
// Function check_parse_result()
// ----------------------------------------------------------------------------
bool check_parse_result(options & options)
{
    if (options.seed_start >= options.seed_end)
    {
        _V(options, "ERROR: Option ''seed_start'' must be lower than ''seed_end''");
        return false;
    }
    if(!options.path_out_files.empty())
    {
        _V(options, "INFO: The specified output file path is " << options.path_out_files);
    }
    else
    {
        std::string tmp_dir;
        tmp_dir = options.tmp_dir;
        if (tmp_dir.empty())
        {
            tmp_dir = SEQAN_TEMP_FILENAME();
            // remove "/test_file" suffix
            tmp_dir.erase(tmp_dir.length() - 10u, tmp_dir.length());
        }
#ifdef STDLIB_VS
        _putenv_s("TMPDIR", tmp_dir.c_str());
#else
        setenv("TMPDIR", tmp_dir.c_str(), true);
#endif
        options.tmp_dir = tmp_dir;
        _V(options, "INFO: The absolute path where to create the tmp_dir is " << tmp_dir);
    }
    if(!options.in_file_primir.empty())
    {
        if(options.in_file_gff.empty())
        {
            _V(options, "ERROR: If primir file is used, a gff file must also be provided");
            return false;
        }
    }
    else if (options.in_file_mature.empty() && options.path_load_serialized.empty())
    {
        _V(options, "ERROR: If primir file is not used, (at least) a mature or a serialized file must be provided");
        return false;
    }

    return true;
}

#endif //ISOMIR_SEA_OPTIONS_H
