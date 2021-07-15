//
//  chromosome.h
//  Generic containers for genetic algorithm
//  Elements that make up chromosomes ("genes") are
//  actually stored elsewhere (in a "gene pool") and
//  "chromosomes" contain justthe index referring
//  to these genes. The methods are therefore just
//  permutations and references to new indexes.
//
//  "Population" class determines a hierarchy between
//  chromosomes to determine which are optimal in
//  reproduction and which are "way better than average"
//
//  Created by Andrea Caprotti on 27/05/21.
//

#ifndef chromosome_h
#define chromosome_h

#include <vector>
#include <algorithm>                // std::sort
#include <iomanip>                  // std::setw
#include <cmath>
#include "random.h"

class Chromosome {
private:
    int gene_no;                    // number of genes
                                    // (dim of vector)
    int ind_min, ind_max;           // saved to remember where chromosome is "cropped"
    int ind, other_ind;             // useful cycle index
    
    std::vector<int> cropped_gene;  // placeholder
    bool can_insert;                // insert only after cropped
    std::vector<int> insertion;

public:
    std::vector<int> gene_sequence; // recognises genes by index
    double weight;                  // fitness index
    
    Chromosome (int N){             // constructor
        gene_no = N;
        
        ind_min = 0;
        ind_max = gene_no;
        
        gene_sequence.clear();      // initialises vector
        cropped_gene.clear();

        can_insert = false;
        
        return;
    }
    
    ~Chromosome(){}                 // destructor
    
    void FillChromo(std::vector<int> &input){ // rapid way of insertion
        gene_sequence.clear();
        for (ind = 0; ind < gene_no; ++ind)
            gene_sequence.push_back(input[ind]);
        return;
        
    }
    
    void PrintSequence(std::ofstream &Out){
        for (ind = 0; ind < gene_no; ind++)
            Out << gene_sequence[ind] << std::setw(3);
        Out << std::endl;
    }
    
    bool GoodChromo (){             // to check if sequence is acceptable
        bool goodness = false;
        double target, sum;
        target = gene_no * (gene_no - 1)/ 2; // expected sum of gene_no-1 int
        sum = 0;
        
        for (ind = 0; ind < gene_no; ++ind)
            sum += gene_sequence[ind];
        if (gene_sequence[0] == 0)           // first gene fixed
            if (sum == target)               // all genes appear exactly once
                goodness = true;
        
        return goodness;
    }
    
    void InsertSection(std::vector<int> &section){
        // inserts new sequence in a previously cropped section
        // inserts only values that have been previosly cropped
        // only works if it has been previously cropped and with right size
        // (checks to keep vector size always constant)
        if (can_insert && ind_max - ind_min == section.size()){
            insertion.clear();              // always initialise!
            // to make sure no indexes are repeated
            
            for ( ind =0; ind <section.size(); ++ind){
                for (other_ind=0; other_ind < cropped_gene.size(); ++other_ind){
                    
                    if (section[ind] == cropped_gene [other_ind]){
                    insertion.push_back(section[ind]);
                    std::swap(cropped_gene[other_ind],cropped_gene.back());
                    cropped_gene.pop_back();
                    // elements are deleted from stored vector to speed-
                    // up following comparisons
                    }
                }
            }
            insertion.insert(insertion.end(), cropped_gene.begin(), cropped_gene.end()); // concatenates new sequence
            
            for (ind = 0; ind < ind_max - ind_min; ++ind) // inserts new seq
                gene_sequence[ind+ind_min] = insertion[ind];
            can_insert = false;
        }
        return;
    }
    
    std::vector<int> CropChro (int i_min, int i_max){
        ind_min = i_min;        // cropping zone assinged by user
        ind_max = i_max;
        
        cropped_gene.clear();
        for (ind = ind_min; ind  < ind_max; ++ind )
            cropped_gene.push_back(gene_sequence [ind]);
        
        can_insert = true;
        return cropped_gene;
        // cropped_gene also saves indexes of chromosome
    }
    
    void SinglePermutation (int ind1, int ind2){
        // swaps genes
        if (ind1 == ind2)   // just to avoid swapping same element…
            ind2 -=1;
        
        if (ind2>0)         // doesn't swap first element
            std::swap (gene_sequence[ind1], gene_sequence[ind2]);
            
        return;
    }
    
    void BlockPermutation(int n_move, int start_ind){
        // receives how many genes are to be moved and starting point
        if (n_move > gene_no/2)
            // otherwise sequence can't be fully moved
            return;
        if (start_ind == 0)
            return;         // to be sure that the first element is fixed
        
        other_ind = start_ind;      // saves starting ind
        while (other_ind > gene_no/2 - n_move)
            other_ind-=1;
        // in order to ensure that the two blocks to be switched surely don't overlap
        
        for (ind = 0; ind < n_move; ++ind)
            std::swap(gene_sequence[other_ind + ind], gene_sequence[other_ind + ind + gene_no/2]);        
        return;
    }
    
    void PartialReverse(int start_ind, int end_ind){
        // reverses genes between given indexes
        std::reverse(gene_sequence.begin() + start_ind, gene_sequence.begin() + end_ind);
    }
    
};

// Population class that collects a bunch of Chromosomes and orders them hierarchically

class PopulationChrome : public Random{
private:
    double skip_child, tot_children, failed_mutations;  // indicators
    double pop_size, gene_no;
    double exponent;
    int ind, other_ind;                                 // cycle indexes
public:
    std::vector<Chromosome> population; // determines hierarchy of
                                        // chromosomes (identified by weight)
    std::vector<Chromosome> store_pop;  // storage for new generation
    
    PopulationChrome(int N, double exp){               // constructor
        gene_no = N;
        skip_child = 0;
        tot_children = 0;
        failed_mutations = 0;
        exponent = exp;
        
        ResetPop();
        return;
    }
    
    ~PopulationChrome(){                                // destructor
        std::cout << skip_child << " skipped children over " << tot_children << " total crossovers " << std::endl;
        std::cout << failed_mutations << " failed attempts of mutation for not following requests" << std::endl;
    }
    
    void ResetPop(){
        population.clear();
        pop_size = 0;
        return;
    }
    
    void AddChromosome (Chromosome chro){
    // best to keep population vector private member,
    // so that ordering is always internal
    // Chromosome needs to be already initialised !!!
        population.push_back(chro);
        pop_size++;
    }
    
    void ReplaceChromo (Chromosome chro, int indexx){
        population[indexx] = chro;
        return;
    }
    
    Chromosome GetBest(){             // returns best element
        return population[0];
    }
    
    double AvgWeight(){        // returns average of best half of population
        double sum = 0;        // as index of total fitness
        for (ind = 0; ind < pop_size; ++ind)
            sum += population[ind].weight;
        return sum/(double)pop_size;
    }
    
    void Sort(){
        // sorts vector hierchically based on weight of each element
        // best has lower weight
        std::sort(population.begin(), population.end(),
                  [](const Chromosome& i, const Chromosome& j) { return i.weight < j.weight; } );
    }
    
    void RandomMutation(int no_mutations, double p0, double p1, double p2){
    // Receives weights for the three possible mutations
    // Monte Carlo choice for acceptance of mutations
        double prob;
        int choice, mut_ind_1, mut_ind_2;
        std::vector<int> placeholder;
        bool goodness;

        for (ind = 0; ind < no_mutations; ++ind){
            other_ind = (int) pop_size*Rannyu();
            // random choice of which chromosome to mutate; all mutate with same prob
            prob = Rannyu();    // acceptance of mutation
            choice = (int)Rannyu(0,3);  // choice of mutation
            
            mut_ind_1 = (int)((gene_no - 1)*Rannyu()) + 1;
            mut_ind_2 = (int)((gene_no - 1)*Rannyu()) + 1;

            // indexes of mutation (avoid first element)
            // second index always greater than first
            
            placeholder = population[other_ind].gene_sequence; //failsafe
            switch (choice) {
                case 0:
                    if (prob < p0)
                        population[other_ind].SinglePermutation (mut_ind_1, mut_ind_2);
                    break;
                    
                case 1:
                    mut_ind_2 = (int)(gene_no/2 * Rannyu());//max genes swapped
                    if (prob < p1)
                        population[other_ind].BlockPermutation (mut_ind_2, mut_ind_1);
                    break;
                    
                default:
                    if (prob < p2){
                        if (mut_ind_1 < mut_ind_2)
                            population[other_ind].PartialReverse (mut_ind_1, mut_ind_2);
                        else
                            population[other_ind].PartialReverse (mut_ind_2, mut_ind_1);
                    }
                    break;
            }
            goodness = population[other_ind].GoodChromo();
            if (!goodness){
                failed_mutations++;
                population[other_ind].FillChromo(placeholder);
            }
            // if mutation doesn't follow requirements it switches back
        }
        return;
    }
    
    void NextGen(double acc_prob){
        // to generate next generation of chromosomes,
        // crossovers are generated and added to storage vector
        // to be replaced afterwards
        
        double prob;        // used to accept crossover
        int ind1, ind2;     // determines section to crossover
        std::vector<int> placeholder1, placeholder2;
        Chromosome first_child(gene_no);
        Chromosome second_child(gene_no);
        bool goodness;
        
        Sort();            // never a bad idea to be sure it's sorted…
        store_pop.clear();
        
        for (int i = 0; i < pop_size/2; ++i){         // generates exactly
                                                      // pop_size new children
            ind1 = (int)((gene_no-1) * Rannyu()) + 1; // extremes of sequence
            ind2 = (int)((gene_no-1) * Rannyu()) + 1; // to swap (avoids 0)
            
            if (ind2 < ind1)
                std::swap (ind1, ind2);
            // new generation given by swap of sequences of genes between two randomly chosen chromosomes with linear probability
            
            ind = (int)(pop_size * pow(Rannyu(),exponent));
            //ind = (int)(pop_size * (1 - sqrt(1 - Rannyu())));
            other_ind = ind;
            while (other_ind == ind) // check for the unlucky time it returns the same index…
                //other_ind = (int)(pop_size * (1 - sqrt(1 - Rannyu())));
                other_ind = (int)(pop_size * pow(Rannyu(),exponent));
            prob = Rannyu();
            
            if (prob < acc_prob){
                first_child = population[ind];
                second_child = population[other_ind];
                
                placeholder1 = first_child.CropChro(ind1,ind2);
                placeholder2 = second_child.CropChro(ind1,ind2);
                
                second_child.InsertSection(placeholder2);
                first_child.InsertSection(placeholder1);
                
                // turns back if offspring doesn't follow requirements
                goodness = first_child.GoodChromo();
                if (!goodness){
                    first_child = population[ind];
                    skip_child++;
                }
                store_pop.push_back(first_child);
                
                goodness = second_child.GoodChromo();
                if (!goodness){
                    second_child = population[other_ind];
                    skip_child++;
                }
                store_pop.push_back(second_child);
                
                tot_children += 2;
            }
            else{                               // crossover only 50% of times
                store_pop.push_back(population[ind]);
                store_pop.push_back(population[other_ind]);
            }
        }
        return;
    }
};

#endif /* chromosome_h */
