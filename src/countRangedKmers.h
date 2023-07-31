#ifndef __COUNTRANGEDKMERS__
#define __COUNTRANGEDKMERS__

std::unordered_map<std::string,int> countRangedKmers(
    std::string sequence,
    std::vector<int> starts,
    std::vector<int> ends,
    int k);

#endif // __COUNTRANGEDKMERS__)