/*
 * =============================================================================
 * iupac_map.h
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_IUPAC_MAP_H
#define ISOMIR_SEA_IUPAC_MAP_H

template<typename t_nt>
using t_key = std::pair<t_nt, t_nt>;

template<typename t_nt>
using t_iupac_map = std::map<t_key<t_nt>, t_nt>;

const t_iupac_map<char> iupac_map =
        {
                {{'A', 'C'}, 'M'},
                {{'A', 'G'}, 'R'},
                {{'A', 'T'}, 'W'}, {{'A', 'U'}, 'W'},
                {{'C', 'G'}, 'S'},
                {{'C', 'T'}, 'Y'}, {{'C', 'U'}, 'Y'},
                {{'G', 'T'}, 'K'}, {{'G', 'U'}, 'K'},
                {{'A', 'S'}, 'V'}, {{'C', 'R'}, 'V'}, {{'G', 'M'}, 'V'},
                {{'A', 'Y'}, 'H'}, {{'C', 'W'}, 'H'}, {{'T', 'M'}, 'H'}, {{'U', 'M'}, 'H'},
                {{'A', 'K'}, 'D'}, {{'G', 'W'}, 'D'}, {{'T', 'R'}, 'D'}, {{'U', 'R'}, 'D'},
                {{'C', 'K'}, 'B'}, {{'G', 'Y'}, 'B'}, {{'T', 'S'}, 'B'}, {{'U', 'S'}, 'B'}
        };

char nt_iupac_converter(char const a, char const b)
{
    std::toupper(a);
    std::toupper(b);
    if(a == b)
    {
        return a;
    } else if (a < b)
    {
        t_key<char> k = {a,b};
        if(iupac_map.count(k)) {
            return iupac_map.at(k);
        }
    } else
    {
        t_key<char> k = {b,a};
        if(iupac_map.count(k)) {
            return iupac_map.at(k);
        }
    }
    return 'N';
}

#endif //ISOMIR_SEA_IUPAC_MAP_H
