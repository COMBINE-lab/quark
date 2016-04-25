// STL
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unordered_map>
#include <tuple>
#include <list>
#include <set>
#include <boost/functional/hash.hpp>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <stack>
#include <queue>
#include <ostream>
// for std::cout

// Boost
typedef std::pair<int,int> Edge;
typedef std::vector<int> Path ;
typedef std::unordered_map< Edge, double, boost::hash< Edge > > EdgeList ;
typedef std::unordered_map< int, std::list<int> > AdjList;
typedef std::set<int> Vertices;
typedef std::unordered_map<int, std::list<std::string> > EqClasses;





int RandomSample(Vertices& vs){
	std::srand(std::time(0));
	int steps = std::rand() % vs.size();
	auto it = vs.begin();
	for(int i =0;i < steps;i++){
		++it;
	}
	int elem = *it;
	//vs.erase(it);
	return elem;
}
//MyTravel(gadj,gedg,vs,randomNode,path);
/**
 DFS(G,v)   ( v is the vertex where the search starts )
         Stack S := {};   ( start with an empty stack )
         for each vertex u, set visited[u] := false;
         push S, v;
         while (S is not empty) do
            u := pop S;
            if (not visited[u]) then
               visited[u] := true;
               for each unvisited neighbour w of u
                  push S, w;
            end if
         end while
      END DFS()
**/
void DFS(AdjList& gadj, EdgeList& gedg, Vertices& vs, Vertices& seenVertices, int randomNode, Path& path_t){
	std::stack<int> dfsstack;
	dfsstack.push(randomNode);
	while(!dfsstack.empty()){
		// print stack for testing
		//std::cout << "\n Stack " <<" \n" ;


		int v=dfsstack.top();
		//std::cout << "curr node :" << v << "\n";
		dfsstack.pop();
		if(seenVertices.find(v) == seenVertices.end()){
			path_t.push_back(v);
			seenVertices.insert(v);
			vs.erase(v);
			std::list<int> neighbors = gadj[v];
			std::vector<Edge> incEdges;
			// take the unvisited neighbors
			// and their edges
			for(auto i: neighbors){
				if(seenVertices.find(i) == seenVertices.end()){
					Edge e = std::make_pair(v,i);
					incEdges.push_back(e);
				}
			}
			// sort the edges according to weight
			std::sort(incEdges.begin(),incEdges.end(),
				[&gedg](const Edge& e1, const Edge& e2) -> bool {
					return gedg[e1] < gedg[e2];
				});
			for(auto &ei : incEdges){
				//std::cout <<"("<< ei.first <<"," << ei.second<<")" <<gedg[ei] << "\n";
				dfsstack.push(ei.second);

			}
		}
	}
}



int main(int argc,char *argv[]){
	if(argc < 2){
	  std::cout << "Give a file name \n" ;
	  return 1;
  	}



  	AdjList gadj;
  	EdgeList gedg;
  	Vertices vs,seenVertices;
  	Path path_t ;

  	std::ifstream is ;
	is.open(argv[1], std::ifstream::in);
    int readCnt = 0;
	for(std::string line; std::getline(is,line); ){
		if(line.empty()) continue;
        if (readCnt % 10000 == 0){
            std::cout << readCnt + 1 << " Edges processed \r" << std::flush ;
        }
        readCnt++;

		int u,v ;
		double w;
		std::istringstream linestream(line);
		linestream >> u >> v >> w;
		gedg[std::make_pair(u,v)] = w;
		gedg[std::make_pair(v,u)] = w;
		//check if the vertex u existed or not
		auto searchu = vs.find(u);
		if(searchu != vs.end()){
			gadj[u].push_back(v);
		}else{
			vs.insert(u);
			gadj[u] = {v};
		}
		//check if the vertex v existed or not
		auto searchv = vs.find(v);
		if(searchv != vs.end()){
			gadj[v].push_back(u);
		}else{
			vs.insert(v);
			gadj[v] = {u};
		}
	}
	int numVertices = vs.size();

	//Graph Loaded
	//We can start doing bfs from random node
	while(seenVertices.size() < numVertices){
		int randomNode = RandomSample(vs);
		//test random sample
		//std::cout<<"random node "<<randomNode<<"\n";
		//sample a random node from vs
		//seenVertices.insert(randomNode);
		//randomNode = 5;
		DFS(gadj,gedg,vs,seenVertices,randomNode,path_t);

		//increase seenVertices;

	}

	EqClasses eqc ;
	std::ifstream ieq ;
	ieq.open(argv[2], std::ifstream::in);
	for(std::string line; std::getline(ieq,line); ){
		int key,numReads;
		std::istringstream linestream(line);
		linestream >> key;
        std::string line2 ;
        std::getline(ieq,line2);
        std::istringstream iss(line2);
        iss >> numReads;
		std::list<std::string> readNames = {};
		while(numReads){
            std::string rname;
            std::getline(ieq,rname);
            readNames.push_back(rname);
            numReads--;
		}
		eqc[key] = readNames ;
	}

	std::ofstream ofs;
	ofs.open(argv[3], std::ofstream::out);
	for(auto i: path_t){
		std::list<std::string> readIds = eqc[i];
		for(auto rid: readIds)
			ofs << rid << "\n";
	}

	/*
	//testing starts
	std::cout << "Iterate through adj list \n";
	for (const auto&u : gadj){
		std::cout<< "[ "<<u.first <<"]:";
		for(const auto&l : u.second)
			std::cout << l << ",";
		std::cout << "\n" ;
	}
	std::cout << "Path \n";
	for (const auto&p : path_t){
		std::cout<<p<<"\n";
	}
	//testing ends
	*/

}
