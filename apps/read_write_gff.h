/*
 * =============================================================================
 * read_write_gff.h
 *
 * Functions to read and write GFF files.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_READ_WRITE_GFF_H
#define ISOMIR_SEA_READ_WRITE_GFF_H

#include <seqan3/range/view/istreambuf.hpp>

inline auto constexpr is_num_sign = is_char<'#'>; // because this is not present in "predicate.hpp" of seqan library
inline auto constexpr is_new_line = is_char<'\n'>; // because this is not present in "predicate.hpp" of seqan library

// ----------------------------------------------------------------------------
// Struct gff_record, adaptation of (old seqan) class GffRecord <include/seqan/gff_io/gff_io_base.h>
// See http://gmod.org/wiki/GFF3 for complete GFF3 file format information
// ----------------------------------------------------------------------------
struct gff_record
{
    //Static member with invalid/sentinel rID value.
    static int32_t const INVALID_POS = 2147483647;

    // The sequence name of the record.
    // The ID of the landmark used to establish the coordinate system for the current feature, most often the
    // contig/chromosome name.
    std::string ref;

    // The source of the record.
    // The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this
    // feature.
    std::string source;

    // The type of the record.
    std::string type;

    // The names of the attributes of the record.
    // For each value there is a name associated in tag_names.
    t_vect_str tag_names;

    // The values of the attributes of the record.
    // For each name there is a value associated in tag_values.
    std::vector<std::string> tag_values;

    // The begin position of the record.
    uint32_t begin_pos;

    // The end position of the record.
    uint32_t end_pos;

    // The score of the record.
    float score;

    // The strand the record belongs to.
    // The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for
    // features that are not stranded.
    char strand;

    // The phase of the record.
    // For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
    // The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the
    // beginning of this feature to reach the first base of the next codon.
    char phase;

    // Returns invalid score (NaN float value).
    // The term <tt>x != x</tt> (for <tt>float x</tt> is only true if <tt>x</tt> is a NaN.
    static float INVALID_SCORE()
    {
        union
        {
            uint32_t u;
            float f;
        } tmp;
        tmp.u = 0x7F800001;
        return tmp.f;
    }
    gff_record() :
            begin_pos(-1), end_pos(-1), score(INVALID_SCORE()),
            strand('.'), phase('.')
    {}
};

// ----------------------------------------------------------------------------
// Function skip_line()
// ----------------------------------------------------------------------------

void skip_line(auto & it, auto & end)
{
    for (; (it != end) && ((!is_new_line)(*it)); ++it) {}
    if(it!= end) { ++it; }
}

// ----------------------------------------------------------------------------
// Function skip_space()
// ----------------------------------------------------------------------------

void skip_space(auto & it, auto & end)
{
    while((it != end) && (is_space(*it)))
    {
        ++it;
    }
}

// ----------------------------------------------------------------------------
// Function skip_empty_lines()
// ----------------------------------------------------------------------------

void skip_empty_lines(auto & it, auto & end)
{
    while((it != end) && (is_space(*it)))
    {
        skip_line(it, end);
    }
}

// ----------------------------------------------------------------------------
// Function skip_comment_directives()
// ----------------------------------------------------------------------------

void skip_comment_directives(auto & it, auto & end)
{
    while((it != end) && (is_num_sign(*it)))
    {
        skip_line(it, end);
    }
}

// ----------------------------------------------------------------------------
// Function read_field()
// ----------------------------------------------------------------------------

void read_field(std::string & field, auto & it, auto & end)
{
    for (; (it != end) && ((!is_new_line)(*it)) && ((!is_blank)(*it)) && *it!='"' && *it!=';'; ++it)
    {
        field.push_back(*it);
    }
    if(it!= end) { ++it; }
}

// ----------------------------------------------------------------------------
// Function read_field_last()
// ----------------------------------------------------------------------------

void read_field_last(std::string & field, auto & it, auto & end)
{
    for (; (it != end) && ((!is_new_line)(*it)) && ((!is_blank)(*it)) && *it!='"' && *it!=';'; ++it)
    {
        field.push_back(*it);
    }
}

// ----------------------------------------------------------------------------
// Function parse_gff_name_value(), adaptation of (old seqan) _parseReadGffKeyValue <include/seqan/gff_io/gff_io_base.h>
// ----------------------------------------------------------------------------

void parse_gff_name_value(std::string & name, std::string & value, auto & it, auto & end)
{
    char c = *it;
    if (c == ' ' || c == '=') {
        debug_stream << "ATTENTION: The key field of an attribute is empty!" << "\n"; // TODO decide how to handle this
    }
    for (; it!=end; ++it)
    {
        c = *it;
        if ((is_new_line)(c) || c == ' ' || c == '=' || c == ';') {
            break;
        }
        name.push_back(c);
    }
    if (it!=end && *it == ';')
    {
        ++it;
        return;
    }
    if ((is_new_line)(*it)) {
        return;
    }
    skip_space(it, end);
    if (*it == '=')
    {
        ++it;
        skip_space(it, end);
    }
    if (*it == '"')
    {
        // Handle the case of a string literal.
        ++it;
        skip_space(it, end);
        read_field(value, it, end);
        // Go over the trailing semicolon and any trailing space.
        for(; it!=end && *it!=';'; ++it) {}
        ++it;
    }
    else
    {
        // Read until the first semicolon, return at whitespace.
        read_field_last(value, it, end);
        if((!is_new_line)(*it)) {
            ++it;
        }
        else {
            return;
        }
        // Skip spaces if any.
        skip_space(it, end);
    }
    return;
}

// ----------------------------------------------------------------------------
// Function read_gff_record(), adaptation of (old seqan) function readRecord() <include/seqan/gff_io/gff_io_base.h>
// ----------------------------------------------------------------------------

void read_gff_record(gff_record & record, auto & it, auto & end)
{
    record = gff_record();
    std::string buffer;
    // Skip empty lines
    skip_empty_lines(it, end);
    // Skip commented line (#) and directives (##)
    skip_comment_directives(it, end);
    // Skip empty lines
    skip_empty_lines(it, end);
    // Read column 1: seqid
    read_field(record.ref, it, end);
    // Read column 2: source
    read_field(record.source, it, end);
    if (record.source == ".") {
        record.source.clear();
    }
    // Read column 3: type
    read_field(record.type, it, end);
    // Read column 4: begin position
    read_field(buffer, it, end);
    record.begin_pos = static_cast<uint32_t>(std::stoul(buffer));
    --record.begin_pos;  // Translate from 1-based to 0-based.
    // Read column 5: end position
    buffer.clear();
    read_field(buffer, it, end);
    record.end_pos = static_cast<uint32_t>(std::stoul(buffer));
    // Check if end < begin
    if (record.end_pos < record.begin_pos) {
        debug_stream << "ATTENTION: Begin position of GFF/GTF record is larger than end position!"
                     << "\n"; // TODO decide how to handle this
    }
    // Read column 6: score
    buffer.clear();
    read_field(buffer, it, end);
    if (buffer != ".") {
        record.score = std::stof(buffer);
    }
    // Read column 7: strand
    buffer.clear();
    read_field(buffer, it, end);
    record.strand = buffer[0];
    // Read column 8: phase
    buffer.clear();
    read_field(buffer, it, end);
    record.phase = buffer[0];
    // It's fine if there are no attributes and the line ends here
    if (it==end || ((is_new_line)(*it)))
    {
        skip_line(it, end);
        return;
    }
    // There is often a space character between phase and attribute columns. Skip that!
    skip_space(it, end);
    // read column 9: attributes
    while (it != end)
    {
        std::string name;
        std::string value;
        // Read next "name=value;" pair.
        parse_gff_name_value(name, value, it, end);
        record.tag_names.push_back(name);
        record.tag_values.push_back(value);
        name.clear();
        value.clear();
        // At end of line: Skip EOL and break.
        if (it!=end && ((is_new_line)(*it)))
        {
            ++it;
            break;
        }
    }
    // The last line might be a "### directive" specifically in GFF3
    // Skip empty lines
    skip_empty_lines(it, end);
    // Skip commented line (#) and directives (##)
    skip_comment_directives(it, end);
    // Skip empty lines
    skip_empty_lines(it, end);
    return;
}

// ----------------------------------------------------------------------------
// Function write_attributes(), adaptation of (old seqan) function _writeAttributes() <include/seqan/gff_io/gff_io_base.h>
// ----------------------------------------------------------------------------

template <typename t_file>
void write_attributes(t_file & out_stream, gff_record const & record) // Works only for GFF file (not GTF)
{
    for (unsigned i = 0; i < record.tag_names.size(); ++i)
    {
        if (i != 0)
        {
            out_stream << ";";
        }
        for(auto it = record.tag_names[i].begin(); it != record.tag_names[i].end(); ++it) {
            if (*it == '\n' || *it == '"') {
                debug_stream << "ATTENTION: Attribute contains illegal character!\n"; // TODO decide how to handle this
            }
            if((*it) == ';' || (*it) == '=')
            {
                out_stream << '"' << record.tag_names[i] << '"';
                break;
            }
        }
        out_stream << record.tag_names[i];

        if (!record.tag_values[i].empty())
        {
            out_stream << "=";
            for(auto it = record.tag_values[i].begin(); it != record.tag_values[i].end(); ++it) {
                if (*it =='\n' || *it == '"') {
                    debug_stream << "ATTENTION: Attribute contains illegal character!\n"; // TODO decide how to handle this
                }
            }
            out_stream << record.tag_values[i];
        }
    }
    return;
}

// ----------------------------------------------------------------------------
// Function write_gff_record(), adaptation of (old seqan) function writeRecord() <include/seqan/gff_io/gff_io_base.h>
// ----------------------------------------------------------------------------

template <typename t_file>
void write_gff_record(t_file & out_stream, gff_record const & record)
{
    // Ignore empty annotations, i.e. annotations that are 'guessed' by implicit information from their children (in GFF)
    //if (record.ref.empty())
    //    return;
    // Write column 1: seqid
    out_stream << record.ref << "\t";
    // Write column 2: source
    if (record.source.empty()) {
        out_stream << "." << "\t";
    }
    else {
        out_stream << record.source << "\t";
    }
    // Write column 3: type
    out_stream << record.type << "\t";
    // Write column 4: begin position
    if (record.begin_pos != (unsigned)-1) {
        out_stream << record.begin_pos + 1 << "\t";
    }
    else {
        debug_stream << "ATTENTION: No start position\n"; // TODO decide how to handle this
    }
    // Write column 5: end position
    if (record.end_pos != (unsigned)-1 && record.begin_pos <= record.end_pos) {
        out_stream << record.end_pos << "\t";
    }
    else {
        debug_stream << "ATTENTION: No end position or begin position larger than end position!\n"; // TODO decide how to handle this
    }
    // Write column 6: score
    if (record.score != record.score) {
        out_stream << "." << "\t";
    }
    else {
        out_stream << record.score << "\t";
    }
    // Write column 7: strand
    out_stream << record.strand << "\t";
    // Write column 8: phase
    out_stream << record.phase << "\t";
    // Write column 9: attributes (only until length - 1, because there is no semicolon at the end of the line)
    write_attributes(out_stream, record);
    out_stream << "\n";
    return;
}

#endif //ISOMIR_SEA_READ_WRITE_GFF_H
