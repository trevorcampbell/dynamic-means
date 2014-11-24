#ifndef __RESULTS_HPP
namespace dmeans{
template<class P>
class Results{
	public:
		std::map<uint64_t, uint64_t> lbls;
		std::map<uint64_t, P> prms;
		double tTaken;
		double cost;
};
}
#define __RESULTS_HPP
#endif /* __RESULTS_HPP */
