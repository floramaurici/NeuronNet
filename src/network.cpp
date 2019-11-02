#include "network.h"
#include "random.h"

void Network::resize(const size_t &n, double inhib) {
    size_t old = size();
    neurons.resize(n);
    if (n <= old) return;
    size_t nfs(inhib*(n-old)+.5);
    set_default_params({{"FS", nfs}}, old);
}

void Network::set_default_params(const std::map<std::string, size_t> &types,
                                 const size_t start) {
    size_t k(0), ssize(size()-start), kmax(0);
    std::vector<double> noise(ssize);
    _RNG->uniform_double(noise);
    for (auto I : types) 
        if (Neuron::type_exists(I.first)) 
            for (kmax+=I.second; k<kmax && k<ssize; k++) 
                neurons[start+k].set_default_params(I.first, noise[k]);
    for (; k<ssize; k++) neurons[start+k].set_default_params("RS", noise[k]);
}

void Network::set_types_params(const std::vector<std::string> &_types,
                               const std::vector<NeuronParams> &_par,
                               const size_t start) {
    for (size_t k=0; k<_par.size(); k++) {
        neurons[start+k].set_type(_types[k]);
        neurons[start+k].set_params(_par[k]);
    }
}

void Network::set_values(const std::vector<double> &_poten, const size_t start) {
    for (size_t k=0; k<_poten.size(); k++) 
        neurons[start+k].potential(_poten[k]);
}

bool Network::add_link(const size_t &a, const size_t &b, double str) {
    if (a==b || a>=size() || b>=size() || str<1e-6) return false;
    if (links.count({a,b})) return false;
    if (neurons[b].is_inhibitory()) str *= -2.0;
    links.insert({{a,b}, str});
    return true;
}

size_t Network::random_connect(const double &mean_deg, const double &mean_streng) {
    links.clear();
    std::vector<int> degrees(size());
    _RNG->poisson(degrees, mean_deg);
    size_t num_links = 0;
    std::vector<size_t> nodeidx(size());
    std::iota(nodeidx.begin(), nodeidx.end(), 0);
    for (size_t node=0; node<size(); node++) {
        _RNG->shuffle(nodeidx);
        std::vector<double> strength(degrees[node]);
        _RNG->uniform_double(strength, 1e-6, 2*mean_streng);
        int nl = 0;
        for (size_t nn=0; nn<size() && nl<degrees[node]; nn++)
            if (add_link(node, nodeidx[nn], strength[nl])) nl++;
        num_links += nl;
    }
    return num_links;
}


std::pair<size_t, double> Network::degree(const size_t& c) const {
	std::vector<std::pair<size_t, double>> connected;
	connected = neighbors(c);
	double sum_intensity(0);
	
	for (size_t i(0); i<connected.size() ; ++i)
		{ sum_intensity += connected[i].second; }
				
	return std::pair<size_t, double> (connected.size(), sum_intensity);	
	}

std::set<size_t> Network::step(const std::vector<double>& tab) {
	std::set<size_t> index;
	double  inhib_sum(0);
	double  exited_sum(0);
	double input_thalamic(0);
	
	for (size_t i(0); i<neurons.size() ; ++i) {
		
		if (neurons[i].firing() ) { 
		  index.insert(i); 
		  neurons[i].reset();
		}
		
		std::vector<std::pair<size_t, double> > connected;
		connected = neighbors(i);
			for (size_t j(0); j<connected.size() ; ++j) {
				
				if (index.count(j) == 1) { 
				if (neurons[j].is_inhibitory() ) { inhib_sum += connected[j].second; }
				else { exited_sum += connected[j].second; }
				}
			}
		
		if (neurons[i].is_inhibitory() ) { input_thalamic = tab[i]*0.4;} //we change the variance when it's an inhibitory neuron
		else {input_thalamic =tab[i]; }
			
		
	
	neurons[i].input( input_thalamic + 0.5*exited_sum + inhib_sum); 
	neurons[i].step();
	
	}
return index;
}  

std::vector<std::pair<size_t, double> > Network::neighbors(const size_t& c) const  
{
	std::vector<std::pair<size_t, double> > neighbor;
	for (std::map<std::pair<size_t, size_t>, double>::const_iterator it=links.lower_bound({c,0}); it!=links.cend() and (it->first.first == c); ++it) { //we begin the iteration at {c,0} to avoid a too long iteration
		
		  	
			std::pair<size_t, double> index (it->first.second, it->second);
			neighbor.push_back(index);
			 
	}
	return neighbor;
}


std::vector<double> Network::potentials() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].potential());
    return vals;
}

std::vector<double> Network::recoveries() const {
    std::vector<double> vals;
    for (size_t nn=0; nn<size(); nn++)
        vals.push_back(neurons[nn].recovery());
    return vals;
}

void Network::print_params(std::ostream *_out) {
    (*_out) << "Type\ta\tb\tc\td\tInhibitory\tdegree\tvalence" << std::endl;
    for (size_t nn=0; nn<size(); nn++) {
        std::pair<size_t, double> dI = degree(nn);
        (*_out) << neurons[nn].formatted_params() 
                << '\t' << dI.first << '\t' << dI.second
                << std::endl;
    }
}

void Network::print_head(const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons)
            if (In.is_type(It.first)) {
                (*_out) << '\t' << It.first << ".v"
                        << '\t' << It.first << ".u"
                        << '\t' << It.first << ".I";
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << "RS.v" << '\t' << "RS.u" << '\t' << "RS.I";
                break;
            }
    (*_out) << std::endl;
}

void Network::print_traj(const int time, const std::map<std::string, size_t> &_nt, 
                         std::ostream *_out) {
    (*_out)  << time;
    size_t total = 0;
    for (auto It : _nt) {
        total += It.second;
        for (auto In : neurons) 
            if (In.is_type(It.first)) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    }
    if (total<size())
        for (auto In : neurons) 
            if (In.is_type("RS")) {
                (*_out) << '\t' << In.formatted_values();
                break;
            }
    (*_out) << std::endl;
}
