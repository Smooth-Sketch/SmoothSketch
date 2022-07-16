#include "basic_sketch/basic_sketch.h"

#include <assert.h>

#ifndef _SMOOTH_BASELINE_H_
#define _SMOOTH_BASELINE_H_


#include <bits/stdc++.h>
#include "Common/hash.h"
#include "Param.h"


template<typename ID_TYPE>
class CountMinSketch{
private:

public:
	uint32_t **CMSketch;
	uint32_t CounterOfArray;

	CountMinSketch(){}

    void init(uint32_t memory){
		CounterOfArray = memory * 1024 / (_NumOfArray * sizeof(uint32_t));
		// std::cout<<"counter of each array: "<<CounterOfArray<<std::endl;
		// CMSketch = new uint32_t*[_NumOfArray];
        CMSketch = (uint32_t**)CALLOC(_NumOfArray, sizeof(uint32_t*));
		for(int i = 0; i < _NumOfArray; ++i){
			// CMSketch[i] = new uint32_t[CounterOfArray]();
            CMSketch[i] = (uint32_t *)CALLOC(CounterOfArray, sizeof(uint32_t));
            memset(CMSketch[i], 0, CounterOfArray * sizeof(uint32_t));
		}
    }
    ~CountMinSketch(){}

	void destroy(){
		for(int i= 0; i < _NumOfArray; ++i){
			FREE(CMSketch[i]);
		}
		FREE(CMSketch);
	}

	int insert(ID_TYPE ItemID){
		for(int i = 0; i < _NumOfArray; ++i){
			uint32_t HashValue = hash_s(ItemID,i);
			// uint32_t HashValue = MurmurHash32((const void*)ItemID.data(),KEY_LEN, 50 + i);
			uint32_t loc = HashValue % CounterOfArray;
			// std::cout<<loc<<","<<CounterOfArray<<std::endl;
			CMSketch[i][loc] ++;
		}
		return 1;
	}

	uint32_t query(ID_TYPE ItemID){
		uint32_t freq_min = UINT32_MAX;
		for(int i = 0; i < _NumOfArray; ++i){
			uint32_t loc = hash_s(ItemID, i) % CounterOfArray;
			// uint32_t loc = MurmurHash32((const void*)ItemID.data(), KEY_LEN, 50 + i) % CounterOfArray;
			uint32_t freq_loc = CMSketch[i][loc];
			freq_min = std::min(freq_min,freq_loc);
		}
		return freq_min;
	}

	void clear(){
		for(int i = 0; i < _NumOfArray; ++i){
			memset(CMSketch[i], 0, CounterOfArray * sizeof(uint32_t));
		}
	}

};

typedef long long DATA_TYPE;

// template<typename ID_TYPE, typename RESULT_ID>
class Baseline: public basic_sketch{
private:
	
	std::set<DATA_TYPE> Smooth_Queue = {};
	uint32_t last_timestamp;
	uint32_t win_cnt;
	uint32_t ArrayCounters;
	uint32_t Cur_timestamp;
	std::vector<std::pair<DATA_TYPE, uint32_t>> result;
    CountMinSketch<DATA_TYPE>* Base[P];

public:

    using basic_sketch::operator new;
    using basic_sketch::operator new[];
    using basic_sketch::operator delete;
    using basic_sketch::operator delete[];
    Baseline(){}

	Baseline(int argc, basic_sketch_string *argv):win_cnt(0), last_timestamp(0){
        init_matrix();

        int memory = argv[0].to_int();
		int MemPerSketch = memory / P;

		ArrayCounters = MemPerSketch * 1024 / (_NumOfArray * sizeof(uint32_t));

		for(int i = 0; i < P; ++i){
		 	// Basic[i] = new CountMinSketch<ID_TYPE>(MemPerSketch);             
            Base[i] = (CountMinSketch<DATA_TYPE>*)CALLOC(1, sizeof(CountMinSketch<DATA_TYPE>));
            Base[i]->init(MemPerSketch);
        }

	}

	~Baseline(){
		for(int i = 0; i < P; ++i){
            Base[i]->destroy();
			FREE(Base[i]);
		}

	}

    Baseline(const basic_sketch_string &s){
		init_matrix();

		size_t tmp = 0;
		const char *ss = s.c_str();

		memcpy(&last_timestamp, ss + tmp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(&win_cnt, ss + tmp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(&ArrayCounters, ss + tmp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(&Cur_timestamp, ss + tmp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		for(int i = 0; i < S; ++i){
			Base[i]->CounterOfArray = ArrayCounters;
			for(int j = 0; j < _NumOfArray; ++j){
				memcpy(&Base[i]->CMSketch[j], ss + tmp, ArrayCounters * sizeof(uint32_t));
				tmp += ArrayCounters * sizeof(uint32_t);
			}
		}
	}

    basic_sketch_string *to_String(){
		char *s1 = (char*)CALLOC(4 + ArrayCounters * _NumOfArray * P, sizeof(uint32_t));
		size_t tmp = 0;

		memcpy(s1 + tmp, &last_timestamp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(s1 + tmp, &win_cnt, sizeof(uint32_t));
		tmp += sizeof(uint32_t);
	
		memcpy(s1 + tmp, &ArrayCounters, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		memcpy(s1 + tmp, &Cur_timestamp, sizeof(uint32_t));
		tmp += sizeof(uint32_t);

		for(int i = 0; i < P; ++i){
			for(int j = 0; j < _NumOfArray; ++j){
				memcpy(s1 + tmp, &Base[i]->CMSketch[j], ArrayCounters * sizeof(uint32_t));
				tmp += ArrayCounters * sizeof(uint32_t);
			}
		}

		basic_sketch_string *bs = new basic_sketch_string(s1, tmp);

		return bs;

    }

	// void insert(ID_TYPE id, uint32_t timestamp){
	// 	if(last_timestamp + window_size < timestamp){//时间切换
	// 		transition();
	// 	}
	// 	Base[win_cnt % P]->insert(id);
	// 	Smooth_Queue.insert(id);
	// }

	basic_sketch_reply *insert(const int &argc, const basic_sketch_string *argv){
		basic_sketch_reply *insert_out = new basic_sketch_reply;

		for(int i = 0; i < argc; ++i){
			basic_sketch_string itemid = argv[i];
			long long id = itemid.to_long_long();
			Cur_timestamp ++;

			if(last_timestamp + window_size < Cur_timestamp){
				transition();
				for(auto it = result.begin(); it != result.end(); it ++){
					insert_out->push_back(it->first);
					insert_out->push_back((long long)it->second);
				}
				result.clear();
			}
			Base[win_cnt % P]->insert(id);
			Smooth_Queue.insert(id);

		}
		
		return insert_out;
	}

	void traverse_query(){
		uint32_t freq[P] = {};
		for(auto i = Smooth_Queue.begin(); i != Smooth_Queue.end(); i++){
			double y[P] = {}, b[K + 1] = {}, z[P] = {};
			int flag = 0;
			DATA_TYPE ItemID = *i;
			for(int j = P - 1, k = win_cnt; j >= 0; j--, k--){
				// freq[j] = Basic[j]->query(&ItemID);
				y[j] = Base[(P + k) % P]->query(ItemID);
				z[j] = y[j];
				if(y[j] == 0){
					flag = 1;
				}
			}

			if(flag){
				continue;
			}
			linear_regressing(z,b);


			if(abs(b[K]) < var_thres){
				continue;
			}

			double error = 0;
			for(int j = 0; j < P; ++j){
				z[j] = y[j];
				for(int k = 0; k <= K; ++k){
					z[j] -= pow(j,k) * b[k];
				}
				error += z[j] * z[j];
			}

			if (error / P <= error_thres){
				// result.emplace_back(std::make_pair(ItemID, win_cnt));
				result.emplace_back(std::make_pair(ItemID,win_cnt));

			}
		}
	}

	void transition(){
		traverse_query();
		win_cnt ++;
		Base[win_cnt % P]->clear();
		Smooth_Queue.clear();
		last_timestamp += window_size;
        // result.clear();
	}

	basic_sketch_reply *flush(){
		traverse_query();
		basic_sketch_reply *flush_out = new basic_sketch_reply;
		for(auto it = result.begin(); it != result.end(); it ++){
			flush_out->push_back(it->first);
			flush_out->push_back((long long)it->second);
		}
		result.clear();

		return flush_out;
	}

	static basic_sketch_reply *Insert(void *o, const int &argc, const basic_sketch_string *argv){
		return ((Baseline *)o)->insert(argc, argv);
	}


	static int command_num(){return 2;}

	static basic_sketch_string command_name(int index){
		basic_sketch_string tmp[] = {"insert", "flush"};
		return tmp[index];
	}

	static basic_sketch_reply *Flush(void *o, const int &argc, const basic_sketch_string *argv){
		return ((Baseline *)o)->flush();
	}

	static basic_sketch_func command(int index){
		basic_sketch_func tmp[] = {(Baseline::Insert), (Baseline::Flush)};
		return tmp[index];
	}

	static basic_sketch_string class_name(){return "Baseline";}

	static int command_type(int index){
		int tmp[] = {0,1};
		return tmp[index];
	}

	static char* type_name(){return "SMOOTH_BS";}
};




#endif