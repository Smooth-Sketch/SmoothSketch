#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#include <bits/stdc++.h>
#include <hash.h>
#include <Mmap.h>
#include "SmoothSketch.h"
#include "CorrectDetector.h"
#include "Strawman.h"
#include "SmoothCU.h"



// #define memory 250

struct CAIDA_Tuple {
    uint64_t timestamp;
    uint64_t id;
};



// typedef uint64_t ID_TYPE;
// typedef uint32_t ID_TYPE;
typedef long long ID_TYPE;
typedef std::multimap<ID_TYPE,std::pair<int, std::pair<uint32_t, uint32_t>>>  ARE_TYPE;

class CAIDABenchmark {
public:
	CAIDABenchmark(std::string PATH, std::string RESFILE) {
        load_result = Load(PATH.c_str());
        dataset = (CAIDA_Tuple*)load_result.start;
        length = load_result.length / sizeof(CAIDA_Tuple);
        // std::cout<<length<<std::endl;

        oFile.open(RESFILE, std::ios::app);
        // assert(!oFile);
        if(!oFile)
            std::cout<<"open file failed"<<std::endl;
        oFile <<"Memory,"<<"ratio,"<<"K,"<<"CellNumber,"<<"WindowsRecord,"<<"a_k_Th,"<<"MSE_Th,"<<"potential_thres,"
        // <<"Correct(B),"<<"Report(B),"<<"GroundTruth(B),"<<"PR(B),"<<"RR(B),"<<"F1(B),"<<"Throughput(B),"
        <<"Correct(S),"<<"Report(S),"<<"GroundTruth(S),"<<"PR(S),"<<"RR(S),"<<"F1(S),"<<"ARE(S),"<<"Throughput(S),"
        // <<"Correct(U),"<<"Report(U),"<<"GroundTruth(U),"<<"PR(U),"<<"RR(U),"<<"F1(U),"<<"Throughput(U),"
        // <<"Correct(CM),"<<"Report(CM),"<<"GroundTruth(CM),"<<"PR(CM),"<<"RR(CM),"<<"F1(CM),"<<"Throughput(CM),"
        // <<"Correct(CF),"<<"Report(CF),"<<"GroundTruth(CF),"<<"PR(CF),"<<"RR(CF),"<<"F1(CF),"<<"Throughput(CF),"
        // <<"Correct(LF),"<<"Report(LF),"<<"GroundTruth(LF),"<<"PR(LF),"<<"RR(LF),"<<"F1(LF),"<<"Throughput(LF),"
        <<std::endl;
    }
	~CAIDABenchmark() {
    }

    void Compare(std::vector<std::pair<ID_TYPE, uint32_t>> &result, 
                 std::vector<std::pair<ID_TYPE, uint32_t>> &truth) {
        // assert(result.size() == 0);
        std::multimap<ID_TYPE, uint32_t> gt;
		for (auto &record : truth) {
			gt.insert(std::make_pair(record.first, record.second));
		}
        
        uint32_t correct = 0, total = truth.size(), all = result.size();
        for (auto &r : result) {
            auto begin = gt.lower_bound(r.first);
			auto end = gt.upper_bound(r.first);
            while (begin != end) {
                if (r.second == begin->second) 
                    correct++;
                begin++;
            }
        }
        if(correct == 0)
            std::cout<<"correct is 0"<<std::endl;

        std::cout << "correct: " << correct << " find: " << all << " truth: " << total << std::endl;
        std::cout << "Precision: " << std::fixed << std::setprecision(4) << 100.0 * correct / all << std::endl;
        std::cout << "Recall:    " << std::fixed << std::setprecision(4) << 100.0 * correct / total << std::endl;
        std::cout << "F1:        " << std::fixed << std::setprecision(4) << 2.0 * correct / (all + total) << std::endl;

        oFile<<correct<<","<<all<<","<<total<<","<<100.0 * correct / all<<","<<100.0 * correct / total<<","<<2.0 * correct / (all + total)<<",";
        //
        gt.clear();
    }
    void Check(std::vector<Report_Slot<ID_TYPE>> &result, std::vector<Report_Slot<ID_TYPE>> &truth) {
        std::cout << result.size() << " " << truth.size() << "\n";
        std::cout << result[0].id << " " << result[0].start_window << " " << result[0].end_window << "\n";
        std::cout << truth[0].id << " " << truth[0].start_window << " " << truth[0].end_window <<"\n";
    }
#ifdef PREDICT_MODE
    void Predict_Check(std::map<uint32_t, std::map<ID_TYPE, uint32_t>> &predict, 
                       std::map<ID_TYPE, std::map<uint32_t, uint32_t>> &history) {
        int correct = 0, total = 0;
        for (auto i : predict) {
            for (auto j : i.second) {
                int predict_result = j.second, truth = history[j.first][i.first];
                if (abs(truth - predict_result) <= 5 || (2 * truth >= predict_result && truth <= 2 * predict_result)) {
                    correct++;
                }
                total++;
            }
        }
        std::cout << "Predict: " << total << " Correct: " << correct << "\n";
        std::cout << "Accuracy: " << std::fixed << std::setprecision(4) << 100.0 * correct / total << "\n";
    }
#endif
    void TopK_Check(std::vector<Report_Slot<ID_TYPE>> &result, std::vector<Report_Slot<ID_TYPE>> &truth){
        std::multimap<ID_TYPE, std::pair<uint32_t, uint32_t>> gt;
        for(auto &record : truth){
            gt.insert(std::make_pair(record.id, std::make_pair(record.start_window, record.end_window)));
        }
        uint32_t correct = 0, total = truth.size(), all = result.size();
        for(auto &res: result){
            auto begin = gt.lower_bound(res.id);
			auto end = gt.upper_bound(res.id);
            while (begin != end) {
                if (res.start_window == begin->second.first && res.end_window == begin->second.second) 
                    correct++;
                begin++;
            } 
        }

        std::cout << "TopK_correct: " << correct << " find: " << all << " truth: " << total << std::endl;
        std::cout << "TopK_Precision: " << std::fixed << std::setprecision(4) << 100.0 * correct / all << std::endl;
        std::cout << "TopK_Recall:    " << std::fixed << std::setprecision(4) << 100.0 * correct / total << std::endl;
        std::cout << "TopK_F1:        " << std::fixed << std::setprecision(4) << 2.0 * correct / (all + total) << std::endl;

        oFile<<correct<<","<<all<<","<<total<<","<<100.0 * correct / all<<","<<100.0 * correct / total<<","<<2.0 * correct / (all + total)<<",";
        gt.clear();
    }

    ARE_TYPE duplicate(std::vector<Report_Slot<ID_TYPE>> &result){
        ARE_TYPE res;
        for(auto &gt_per: result){
            int smooth_len = gt_per.end_window - gt_per.start_window;
            res.insert(std::make_pair(gt_per.id, std::make_pair(smooth_len, std::make_pair(gt_per.start_window, gt_per.end_window))));
        }

        for (auto i = res.begin(); i != res.end(); i++){
        auto count = res.count(i->first)-1;
        auto it = res.find(i->first);
        if(count > 0) it++;
        while (count){
            i->second.first += it->second.first;
            // i->second.second.first += it->second.second.first;
            auto tmp = it;
            it++;
            count--;
            res.erase(tmp);
        }
        }
        std::cout<<"after size: "<<res.size()<<std::endl;
        return res;
    }

    float ARE_cal(std::vector<Report_Slot<ID_TYPE>> &result, std::vector<Report_Slot<ID_TYPE>> &truth){
        std::multimap<ID_TYPE,std::pair<int, std::pair<uint32_t, uint32_t>>> gt;
        std::multimap<ID_TYPE,std::pair<int, std::pair<uint32_t, uint32_t>>> smooth_res;

        gt = duplicate(truth);
        smooth_res = duplicate(result);


        int same_id_count = 0;
        float AAE_2 = 0;
        float ARE_2 = 0;
        for (auto i = smooth_res.begin(); i != smooth_res.end(); i++){
            auto it = gt.find(i->first);
            auto count = gt.count(i->first);
                        
            if(count > 0){
                int AAE_per_2 = abs(i->second.first - it->second.first);
                float ARE_per_2 = abs(i->second.first - it->second.first)/it->second.first;
                same_id_count ++;
                AAE_2 += AAE_per_2;
                ARE_2 += ARE_per_2;
            }                            
        }
        std::cout<<"same id count: "<<same_id_count<<std::endl;
        std::cout<<"AAE: "<<AAE_2/same_id_count<<std::endl;
        float ARE = ARE_2/same_id_count;
        return ARE;

    }

    void Run() {
        uint32_t run_length = 30000000;
        // var_thres = var_input;
        // for(var_thres = 0.25; var_thres <= 2; var_thres += 0.25){

            CorrectDetector<ID_TYPE>* correct_detector = new CorrectDetector<ID_TYPE>(error_thres,var_thres);

            for (int i = 0; i < run_length; ++i) 
                correct_detector->insert(dataset[i%length].id, i);
            
            std::vector<std::pair<ID_TYPE, uint32_t>> ground_truth = correct_detector->query();
            std::vector<Report_Slot<ID_TYPE>> truth_top_k = correct_detector->report();
        #ifdef PREDICT_MODE
            std::map<ID_TYPE, std::map<uint32_t, uint32_t>> history = correct_detector->get_history();
        #endif

            // for(potential_thres = 0; potential_thres <= 2; potential_thres += 0.25){
            for(uint32_t memory = 10; memory <= 10; memory += 10){
                std::cout<<"Mem: "<<memory<<std::endl;

                oFile<<memory<<","<<ratio<<","<<K<<","<<bucket_size<<","<<NumOfWin<<","<<var_thres<<","<<error_thres<<","<<potential_thres<<",";

            /*
            {
                Baseline<ID_TYPE,ID_TYPE>* smooth_sketch = new Baseline<ID_TYPE,ID_TYPE>(memory,var_thres,error_thres);
                auto start = std::chrono::high_resolution_clock::now();

                for (int i = 0; i < run_length; ++i) 
                    smooth_sketch->insert(dataset[i%length].id, i);

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::milli>tm = end - start;
    
                std::vector<std::pair<ID_TYPE, uint32_t>> smooth_result = smooth_sketch->query();
                std::vector<Report_Slot<ID_TYPE>> smooth_report = smooth_sketch->report();
                std::cout<<smooth_result.size()<<std::endl;
                std::cout<<"Baseline:"<<std::endl;
                Compare(smooth_result, ground_truth);
                // Check(smooth_report, truth_top_k);
                // TopK_Check(smooth_report, truth_top_k);

                float ARE = ARE_cal(smooth_report, truth_top_k);
                oFile<<ARE<<",";
                std::cout<<"ARE: "<<ARE<<std::endl;

                oFile<<run_length / (1.0 * tm.count())<<",";
            }
            */
            
            {
                
                SmoothSketch<ID_TYPE>* smooth_sketch = new SmoothSketch<ID_TYPE>(memory, var_thres, error_thres, ratio, bucket_size, NumOfWin, potential_thres);
                auto start = std::chrono::high_resolution_clock::now();

                for (int i = 0; i < run_length; ++i) 
                    smooth_sketch->insert(dataset[i%length].id, i);

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::milli>tm = end - start;
    
                std::vector<std::pair<ID_TYPE, uint32_t>> smooth_result = smooth_sketch->query();
                std::vector<Report_Slot<ID_TYPE>> smooth_report = smooth_sketch->report();
                // std::cout<<smooth_result.size()<<std::endl;
                std::cout<<"SmoothSketch:"<<std::endl;
                Compare(smooth_result, ground_truth);
            #ifdef PREDICT_MODE
                std::map<uint32_t, std::map<ID_TYPE, uint32_t>> predict = smooth_sketch->predict();
                Predict_Check(predict, history);
            #endif
                // Check(smooth_report, truth_top_k);
                // TopK_Check(smooth_report, truth_top_k);

                float ARE = ARE_cal(smooth_report, truth_top_k);
                oFile<<ARE<<",";
                std::cout<<"ARE: "<<ARE<<std::endl;

                oFile<<run_length / (1.0 * tm.count())<<std::endl;
            }
            
            /*
            {
                SmoothSketch_CU<ID_TYPE>* smooth_sketch = new SmoothSketch_CU<ID_TYPE>(memory, var_thres, error_thres, ratio, bucket_size, NumOfWin, potential_thres);
                auto start = std::chrono::high_resolution_clock::now();

                for (int i = 0; i < run_length; ++i) 
                    smooth_sketch->insert(dataset[i%length].id, i);

                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::milli>tm = end - start;
    
                std::vector<std::pair<ID_TYPE, uint32_t>> smooth_result = smooth_sketch->query();
                std::vector<Report_Slot<ID_TYPE>> smooth_report = smooth_sketch->report();
                std::cout<<smooth_result.size()<<std::endl;
                std::cout<<"SmoothSketchCU:"<<std::endl;
                Compare(smooth_result, ground_truth);
                // float ARE = ARE_cal(smooth_report, truth_top_k);
                // oFile<<ARE<<",";
                // std::cout<<"ARE: "<<ARE<<std::endl;
                
                // TopK_Check(smooth_report, truth_top_k);

                oFile<<run_length / (1.0 * tm.count())<<",";
            }
            
            */
            
            
            }
            
            oFile<<std::endl;
            // }
        // }
        oFile.close();
    }
private:
	std::string filename;
    LoadResult load_result;
    CAIDA_Tuple *dataset;
    uint64_t length;
    std::ofstream oFile;

    double var_thres = var_thres_p;
    double error_thres = error_thres_p;
    int bucket_size = bucket_size_p;
    int NumOfWin = S_p;
    double ratio = stage_ratio_p;
    double potential_thres = potential_thres_p;

};

#endif
