#include "art/Persistency/Common/Wrapper.h"

#include "SummaryData/summary.h"

template class std::vector<sumdata::RunData>;
template class art::Wrapper<sumdata::RunData>;
template class std::vector<sumdata::POTSummary>;
template class art::Wrapper<sumdata::POTSummary>;

