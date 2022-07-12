#ifndef _UTILS_H_
#define _UTILS_H_


template<class bidiiter>
bidiiter random_choose(bidiiter begin, bidiiter end, size_t num_random) {
    size_t left = std::distance(begin, end);
    while (num_random--) {
        bidiiter r = begin;
		uniform_int_distribution<int> uniform_int_path(0, left - 1);
		int mov = uniform_int_path(generator);
        std::advance(r, mov);
        std::swap(*begin, *r);
        ++begin;
        --left;
    }
    return begin;
}


#endif

