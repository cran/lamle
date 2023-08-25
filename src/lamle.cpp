// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"
#include "Rmath.h"
#include "cmath"
#include "chrono"
#include "Rcpp.h"
#include "string"

using namespace Rcpp;
using namespace arma;

const double log2pi = std::log(2.0 * M_PI);

//[[Rcpp::export]]
double rcpp_factorial(double n) {
	double fac = 1;
	for(int i = 1; i <= n; ++i) {
		fac = i*fac;
	}
	return(fac); // to compute n!
}

//For filtering.
//[[Rcpp::export]]
int find_row ( arma::mat cand_mat, arma::mat pool_mat ) {
    
    //  This function finds that the row number of the matrix pool_mat that has the same elements as cand_mat.
    //  Negative value means that cand_mat is not found.
    //  Otherwise, the row number in pool_mat is returned.
    int ret = -1 ;
    int row_no = pool_mat.n_rows ;
    int col_no = pool_mat.n_cols ;
    
    for( int r_i=0; r_i<row_no; r_i++ ){

        //  Go through all elements to check whether we find the right row number.
        bool same = true ;
        for( int loop=0; loop<col_no; loop++ ){
            if( cand_mat(0,loop) != pool_mat(r_i,loop) ){
                same = false ;
            }
        }

        //  Check whether cand_mat has been found. 
        if( same == true ){
            ret = r_i ;
            break ;
        }
        
    }
    
    return(ret) ;
    
}

//[[Rcpp::export]]
List FindUniqComb_jlmr_GLLVM( arma::mat loadmat, bool check_zero ){
    
    //  loadmat is a block diagonal matrix of beta.y and beta.x.
    int noK_yx = loadmat.n_cols ; //  Total number of columns
    int noJ_yx = loadmat.n_rows ; //  Total number of columns
    
    //  First, we need to vectorize all unique combinations of (i1,i2,i3,i4) and of (i1,i2) into column matrices.
    mat uniq_comb_2 = zeros<mat>(noK_yx,noK_yx) ;
    mat uniq_comb_2_vech = zeros<mat>(0.5*noK_yx*(noK_yx+1),2) ;
    mat uniq_comb_4 = zeros<mat>(1,5) ;
    int total_row_2 = 0 ;
    int total_row_4 = -1 ;
    for( int i1=0; i1<noK_yx; i1++ ){
        for( int i2=i1; i2<noK_yx; i2++ ){
            Rcpp::checkUserInterrupt();
            uniq_comb_2(i1,i2) = total_row_2 ;
            uniq_comb_2(i2,i1) = total_row_2 ;
            uniq_comb_2_vech(total_row_2,0) = i1 ;
            uniq_comb_2_vech(total_row_2,1) = i2 ;
            total_row_2 ++ ;
            
            for( int i3=i2; i3<noK_yx; i3++ ){
                for( int i4=i3; i4<noK_yx; i4++ ){
                    
                    mat tmp_mat_4 = zeros<mat>(1,5) ;
                    tmp_mat_4(0,0) = i1 ;
                    tmp_mat_4(0,1) = i2 ;
                    tmp_mat_4(0,2) = i3 ;
                    tmp_mat_4(0,3) = i4 ;
                    total_row_4 ++ ;
                    tmp_mat_4(0,4) = total_row_4 ;
                    uniq_comb_4 = join_cols( uniq_comb_4,tmp_mat_4 ) ;	
                    
                }
            }
        }
    }
    //  The first row needs to be removed, because we use the first row to initiate the matrix.
    //    total_row_4 already starts from 0. So we need to have total_row_4+1 here!!!! Otherwise, the combination (noK_yx,noK_yx,noK_yx,noK_yx) will be missing.
    uniq_comb_4 = uniq_comb_4.submat( 1,0,total_row_4+1,4 ) ;
    
    //  The term with the forth order partial derivative is of the form (i1,i2,i3,i4)*(i1,i3)*(i2,i4). So we will find the unique combinations first.
    //    (i1,i2,i3,i4), (i1,i3), and (i2,i4) will be indexed by a single number!
    //    At the same time, note that ONLY (eta,eta,eta,eta) and (xi,xi,xi,xi) contain nonzero elements. So only these potential nonzero combinations will be kept.
    mat uniq_index_4 = ones<mat>(1,3) ;  
    mat freq = ones<mat>(1,1) ;
    double prod_load = 0.0 ;
    for( int i1=0; i1<noK_yx; i1++ ){
        for( int i2=0; i2<noK_yx; i2++ ){
			Rcpp::checkUserInterrupt();
            for( int i3=0; i3<noK_yx; i3++ ){
                for( int i4=0; i4<noK_yx; i4++ ){
                    
                    //  For a CFA model, we can filter out zero derivatives.
                    if( check_zero == true ){
                        
                        for ( int j=0; j<noJ_yx; j++ ) {
                            
                            mat beta_j = trans( loadmat.row(j) ) ;
                            prod_load = beta_j(i1,0)*beta_j(i2,0)*beta_j(i3,0)*beta_j(i4,0) ;
                            
                            if( prod_load != 0.0 ) {
                                break ;
                            }
                            
                        }
                        
                    }
                    
                    if( prod_load==0.0 ){
                        //  No need to keep this combination.
                        continue ;
                    }
                    
                    //  We need to include this combination. So we need to identity the index for the combination (i1,i2,i3,i4) and the related combinations (i1,i3) and (i2,i4).
                    //    Sort (i1,i2,i3,i4) into the matrix ghost_mat1.
                    mat tmp_mat_4 = zeros<mat>(1,4) ;
                    tmp_mat_4(0,0) = i1 ;
                    tmp_mat_4(0,1) = i2 ;
                    tmp_mat_4(0,2) = i3 ;
                    tmp_mat_4(0,3) = i4 ;
                    mat ghost_mat1 = tmp_mat_4 ;
                    std::sort( ghost_mat1.begin(), ghost_mat1.end() ) ;
                    //   Find the position of the sorted (i1,i2,i3,i4) in the index matrix uniq_comb_4
                    int row_no = 0 ;
                    for( int r_i=0; r_i<=total_row_4; r_i++ ){
                        mat target_mat = uniq_comb_4.submat(r_i,0,r_i,3) ;
                        
                        int same_1 = 1 ;
                        for( int loop=0; loop<4; loop++ ){
                            if( ghost_mat1(0,loop)!=target_mat(0,loop) ){
                                same_1 = 0 ;
                            }
                        }
                        //  If same is still 1, then this combination has been found. 
                        if( same_1==1 ){
                            row_no = uniq_comb_4(r_i,4) ;
                            break ;
                        }
                    }
                    mat pos_mat = zeros<mat>(1,3) ;
                    pos_mat(0,0) = row_no ;
                    
                    // (i1,i3) and (i2,i4) can be mapped into the matrix uniq_comb_2, since uniq_comb_2 is a noK_yx*noK_yx matrix. 
                    mat all_comb = tmp_mat_4 ;
                    //  Index of (i1,i3) and (i2,i4), sorted from smallest to largest.
                    if( uniq_comb_2( all_comb(0,0),all_comb(0,2) )<=uniq_comb_2( all_comb(0,1),all_comb(0,3) ) ){
                        pos_mat(0,1) = uniq_comb_2( all_comb(0,0),all_comb(0,2) ) ; //  (i1,i3)
                        pos_mat(0,2) = uniq_comb_2( all_comb(0,1),all_comb(0,3) ) ; //  (i2,i4)
                    }
                    else {
                        pos_mat(0,1) = uniq_comb_2( all_comb(0,1),all_comb(0,3) ) ; //  (i1,i3)
                        pos_mat(0,2) = uniq_comb_2( all_comb(0,0),all_comb(0,2) ) ; //  (i2,i4)
                    }
                    
                    //  Now we need to find the uniuqe combinations and their frequencies.
                    mat ghost_mat = pos_mat ;
                    //    We search the combinations from uniq_index_4 to see whether the provided combination exists or not.
                    int maxsearch = uniq_index_4.n_rows ;
                    int same_2 = 1 ;
                    for( int j=0; j<maxsearch; j++ ){
                        
                        //  Check whether we have computed this.
                        mat target_mat = uniq_index_4.row(j) ;
                        same_2 = 1 ;
                        for( int loop=0; loop<3; loop++ ){
                            
                            if( ghost_mat(0,loop)!=target_mat(0,loop) ){
                                same_2 = 0 ;
                            }
                            
                        }
                        //  This combination already exists.
                        if( same_2==1 ){
                            freq(j,0) = freq(j,0) + 1 ;
                            break ;
                        }
                        
                    }
                    
                    //  If same_2 is still 0, then we have a new combination.
                    if( same_2==0 ){
                        uniq_index_4 = join_cols( uniq_index_4,ghost_mat ) ;
                        freq = join_cols( freq,ones<mat>(1,1) ) ;
                    }
                    
                }
                
            }
            
        }
        
    }
    
    uniq_index_4 = uniq_index_4.submat( 1,0,uniq_index_4.n_rows-1,2 ) ;
    freq = freq.submat( 1,0,freq.n_rows-1,0 ) ;
    
    //  Now we will convert the single-number index back to the combination of the matrix row number of column number.
    //  This is done to the Hessian matrix.
    int uniqpat = uniq_index_4.n_rows ;
    mat index_mat = zeros<mat>(uniqpat,5) ;
    index_mat.col(0) = uniq_index_4.col(0) ;
    for ( int j=0; j<uniqpat; j++ ) {
        index_mat(j,1) = uniq_comb_2_vech( uniq_index_4(j,1), 0 ) ;
        index_mat(j,2) = uniq_comb_2_vech( uniq_index_4(j,1), 1 ) ;
        index_mat(j,3) = uniq_comb_2_vech( uniq_index_4(j,2), 0 ) ;
        index_mat(j,4) = uniq_comb_2_vech( uniq_index_4(j,2), 1 ) ;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("UniqComb") = join_rows(index_mat,freq) ,
        Rcpp::Named("Uniq4") = uniq_comb_4
    ) ;
    
}

//[[Rcpp::export]]
List FindUniqComb_jlmrst_GLLVM( arma::mat loadmat, bool check_zero ){
    
    //  This function returns the unique combinations of partial derivatives in the second-order Laplace approximation.
    //  Two 3rd-order derivative related sums are returned separately.
    //  We can also return them together. But it may even be more combinations. 
    int noK_yx = loadmat.n_cols ; //  Total number of columns
    int noJ_yx = loadmat.n_rows ; //  Total number of columns
    
    //  First, we need to vectorize all unique combinations of (i1,i2,i3,i4) and of (i1,i2) into column matrices.
    mat uniq_comb_2 = zeros<mat>(noK_yx,noK_yx) ;
    mat uniq_comb_2_vech = zeros<mat>(0.5*noK_yx*(noK_yx+1),2) ;
    mat uniq_comb_3 = zeros<mat>(1,4) ;
    int total_row_2 = 0 ;
    int total_row_3 = -1 ;
    for( int i1=0; i1<noK_yx; i1++ ){
		Rcpp::checkUserInterrupt();
        for( int i2=i1; i2<noK_yx; i2++ ){
            
            uniq_comb_2(i1,i2) = total_row_2 ;
            uniq_comb_2(i2,i1) = total_row_2 ;
            uniq_comb_2_vech(total_row_2,0) = i1 ;
            uniq_comb_2_vech(total_row_2,1) = i2 ;
            total_row_2 ++ ;
            
            for( int i3=i2; i3<noK_yx; i3++ ){
                
                mat tmp_mat_3 = zeros<mat>(1,4) ;
                tmp_mat_3(0,0) = i1 ;
                tmp_mat_3(0,1) = i2 ;
                tmp_mat_3(0,2) = i3 ;
                total_row_3 ++ ;
                tmp_mat_3(0,3) = total_row_3 ;
                uniq_comb_3 = join_cols( uniq_comb_3,tmp_mat_3 ) ;	
                
            }
        }
    }
    //  The first row needs to be removed, because we use the first row to initiate the matrix.
    //    total_row_3 already starts from 0. So we need to have total_row_3+1 here!!!! Otherwise, the combination (noK_yx,noK_yx,noK_yx) will be missing.
    uniq_comb_3 = uniq_comb_3.submat( 1,0,total_row_3+1,3 ) ; 

    //  The term with the third order partial derivative is of the form (i1,i2,i3)*(i4,i5,i6)*[ (i1,i2)*(i3,i4)*(i5,i6) + (i1,i4)*(i2,i5)*(i3,i6) ]. So we will find the unique combinations first.
    //    All combinations will be indexed by a single number!
    mat uniq_index_3_4 = -1.0*ones<mat>(1,5) ;  
    mat uniq_index_3_6 = -1.0*ones<mat>(1,5) ;  
    mat freq_4 = ones<mat>(1,1) ;
    mat freq_6 = ones<mat>(1,1) ;
    for( int i1=0; i1<noK_yx; i1++ ){
        for( int i2=0; i2<noK_yx; i2++ ){
            for( int i3=0; i3<noK_yx; i3++ ){
                Rcpp::checkUserInterrupt();
                //  Sort i1,i2,i3 from smallest to largest. So the arrangement is (eta,eta,eta), (eta,xi,xi) or (xi,xi,xi).
                //  Then we don't need to consider (xi,eta,xi) and (xi,xi,eta).
                mat tmp_1 = zeros<mat>(1,3) ;
                tmp_1(0,0) = i1 ;
                tmp_1(0,1) = i2 ;
                tmp_1(0,2) = i3 ;
                mat ghost_mat1 = tmp_1 ;
                std::sort( ghost_mat1.begin(), ghost_mat1.end() ) ;
                int ii1 = ghost_mat1(0,0) ;
                int ii2 = ghost_mat1(0,1) ;
                int ii3 = ghost_mat1(0,2) ;
                
                bool proceed123 = true ;
                if( check_zero == true ){
                    
                    //  We will go through all items from the loading matrix to filter out zero derivatives due to the loading structure.
                    //  The implementation is valid for CFA/graded response model.
                    double prod_load = 0.0 ;
                    for ( int j=0; j<noJ_yx; j++ ) {
                        
                        mat beta_j = trans( loadmat.row(j) ) ;
                        prod_load = beta_j(ii1,0)*beta_j(ii2,0)*beta_j(ii3,0) ;
                        
                        if( prod_load != 0.0 ) {
                            break ;
                        }				
                    }
                    
                    if( prod_load==0.0 ){
                        proceed123 = false ;
                    }
                    
                }
                
                //  If the derivative is zero, we don't need to keep this combination.
                if( proceed123 == false ){
                    continue ;
                }
                
                //  Find the position of the sorted (i1,i2,i3) in the index matrix uniq_comb_3
                int row_no_123 = find_row( ghost_mat1, uniq_comb_3.submat(0,0,total_row_3,2) ) ;
                
                for( int i4=0; i4<noK_yx; i4++ ){
                    for( int i5=0; i5<noK_yx; i5++ ){
                        for( int i6=0; i6<noK_yx; i6++ ){
                            
                            //  Sort i4,i5,i6 from smallest to largest.
                            mat tmp_2 = zeros<mat>(1,3) ;
                            tmp_2(0,0) = i4 ;
                            tmp_2(0,1) = i5 ;
                            tmp_2(0,2) = i6 ;
                            mat ghost_mat2 = tmp_2 ;
                            std::sort( ghost_mat2.begin(), ghost_mat2.end() ) ;
                            int ii4 = ghost_mat2(0,0) ;
                            int ii5 = ghost_mat2(0,1) ;
                            int ii6 = ghost_mat2(0,2) ;
                            
                            bool proceed456 = true ;
                            if( check_zero == true ){
                                
                                //  We will go through all items from the loading matrix to filter out zero derivatives due to the loading structure.
                                //  The implementation is valid for CFA/graded response model.
                                double prod_load = 0.0 ;
                                for ( int j=0; j<noJ_yx; j++ ) {
                                    
                                    mat beta_j = trans( loadmat.row(j) ) ;
                                    prod_load = beta_j(ii4,0)*beta_j(ii5,0)*beta_j(ii6,0) ;
                                    
                                    if( prod_load != 0.0 ) {
                                        break ;
                                    }				
                                }
                                if( prod_load==0.0 ){
                                    proceed456 = false ;
                                }
                                
                            }
                            
                            //  If the derivative is zero, we don't need to keep this combination.
                            if( proceed456 == false ){
                                continue ;
                            }
                            
                            //  Find the position of the sorted (i4,i5,i6) in the index matrix uniq_comb_3
                            int row_no_456 = find_row( ghost_mat2, uniq_comb_3.submat(0,0,total_row_3,2) ) ;
                            
                            //  Sort the single-value index of (i1,i2,i3) and (i4,i5,i6)
                            mat pos_mat = ones<mat>(1,2) ;
                            if( row_no_123<=row_no_456 ){
                                pos_mat(0,0) = row_no_123 ;
                                pos_mat(0,1) = row_no_456 ;
                            }
                            else {
                                pos_mat(0,0) = row_no_456 ;
                                pos_mat(0,1) = row_no_123 ;
                            }
                            
                            // (i1,i2) (i3,i4) (i5,i6)  can be mapped into the matrix uniq_comb_2, since uniq_comb_2 is a noK*noK matrix. 
                            mat all_comb = zeros<mat>(1,6) ;
                            all_comb(0,0) = i1 ;
                            all_comb(0,1) = i2 ;
                            all_comb(0,2) = i3 ;
                            all_comb(0,3) = i4 ;
                            all_comb(0,4) = i5 ;
                            all_comb(0,5) = i6 ;
                            //  Index of (i1,i2) (i3,i4) (i5,i6), sorted from smallest to largest.
                            mat tmp_3 = zeros<mat>(1,3) ;
                            tmp_3(0,0) = uniq_comb_2( all_comb(0,0),all_comb(0,1) ) ;
                            tmp_3(0,1) = uniq_comb_2( all_comb(0,2),all_comb(0,3) ) ;
                            tmp_3(0,2) = uniq_comb_2( all_comb(0,4),all_comb(0,5) ) ; 
                            mat ghost_mat3 = tmp_3 ;
                            std::sort( ghost_mat3.begin(), ghost_mat3.end() ) ;
                            mat pos2_mat = ghost_mat3 ;  //  Positions of (i1,i2) (i3,i4) (i5,i6).
                            
                            //  Now we need to find the uniuqe combinations and their frequencies.
                            mat ghost_mat4 = join_rows(pos_mat,pos2_mat) ;
                            //  We search the combinations from uniq_index_3_4 to see whether the provided combination exists or not.
                            //  If the returned row number is negative, it means that this pattern does not exist in the pool.
                            //  Otherwise, it return the corresponding row number.
                            int index_3_4 = find_row( ghost_mat4, uniq_index_3_4 ) ;
                            if( index_3_4 < 0 ){
                                //  We have a new combination
                                uniq_index_3_4 = join_cols( uniq_index_3_4,ghost_mat4 ) ;
                                freq_4 = join_cols( freq_4,ones<mat>(1,1) ) ;
                            }
                            else {
                                freq_4(index_3_4,0) = freq_4(index_3_4,0) + 1 ;
                            }
                            
                            // (i1,i4) (i2,i5) (i3,i6)  can be mapped into the matrix uniq_comb_2, since uniq_comb_2 is a noK_yx*noK_yx matrix. 
                            //  Index of (i1,i4) (i2,i5) (i3,i6), sorted from smallest to largest.
                            mat tmp_4 = zeros<mat>(1,3) ;
                            tmp_4(0,0) = uniq_comb_2( all_comb(0,0),all_comb(0,3) ) ;
                            tmp_4(0,1) = uniq_comb_2( all_comb(0,1),all_comb(0,4) ) ;
                            tmp_4(0,2) = uniq_comb_2( all_comb(0,2),all_comb(0,5) ) ; 
                            mat ghost_mat5 = tmp_4 ;
                            std::sort( ghost_mat5.begin(), ghost_mat5.end() ) ;
                            mat pos3_mat = ghost_mat5 ;  //  Positions of (i1,i4) (i2,i5) (i3,i6). 
                            //  Now we need to find the uniuqe combinations and their frequencies.
                            mat ghost_mat6 = join_rows(pos_mat,pos3_mat) ;
                            //  We search the combinations from uniq_index_3_6 to see whether the provided combination exists or not.
                            //  If the returned row number is negative, it means that this pattern does not exist in the pool.
                            //  Otherwise, it return the corresponding row number.
                            int index_3_6 = find_row( ghost_mat6, uniq_index_3_6 ) ;
                            if( index_3_6 < 0 ){
                                //  We have a new combination
                                uniq_index_3_6 = join_cols( uniq_index_3_6, ghost_mat6 ) ;
                                freq_6 = join_cols( freq_6,ones<mat>(1,1) ) ;
                            }
                            else {
                                freq_6(index_3_6,0) = freq_6(index_3_6,0) + 1 ;
                            }
                            
                        }
                    }
                }
                
            }
        }
    }
    
    uniq_index_3_4 = uniq_index_3_4.submat( 1,0,uniq_index_3_4.n_rows-1,4 ) ;
    uniq_index_3_6 = uniq_index_3_6.submat( 1,0,uniq_index_3_6.n_rows-1,4 ) ;
    freq_4 = freq_4.submat( 1,0,freq_4.n_rows-1,0 ) ;
    freq_6 = freq_6.submat( 1,0,freq_6.n_rows-1,0 ) ;
    
    //  Now we will convert the single-number index back to the combination of the matrix row number of column number.
    //  This is done to the Hessian matrix.
    int uniqpat_3_4 = uniq_index_3_4.n_rows ;
    mat index_mat_3_4 = zeros<mat>(uniqpat_3_4,8) ;
    index_mat_3_4.col(0) = uniq_index_3_4.col(0) ;
    index_mat_3_4.col(1) = uniq_index_3_4.col(1) ;
    for ( int j=0; j<uniqpat_3_4; j++ ) {
        index_mat_3_4(j,2) = uniq_comb_2_vech( uniq_index_3_4(j,2), 0 ) ;
        index_mat_3_4(j,3) = uniq_comb_2_vech( uniq_index_3_4(j,2), 1 ) ;
        index_mat_3_4(j,4) = uniq_comb_2_vech( uniq_index_3_4(j,3), 0 ) ;
        index_mat_3_4(j,5) = uniq_comb_2_vech( uniq_index_3_4(j,3), 1 ) ;
        index_mat_3_4(j,6) = uniq_comb_2_vech( uniq_index_3_4(j,4), 0 ) ;
        index_mat_3_4(j,7) = uniq_comb_2_vech( uniq_index_3_4(j,4), 1 ) ;
    }
    int uniqpat_3_6 = uniq_index_3_6.n_rows ;
    mat index_mat_3_6 = zeros<mat>(uniqpat_3_6,8) ;
    index_mat_3_6.col(0) = uniq_index_3_6.col(0) ;
    index_mat_3_6.col(1) = uniq_index_3_6.col(1) ;
    for ( int j=0; j<uniqpat_3_6; j++ ) {
        index_mat_3_6(j,2) = uniq_comb_2_vech( uniq_index_3_6(j,2), 0 ) ;
        index_mat_3_6(j,3) = uniq_comb_2_vech( uniq_index_3_6(j,2), 1 ) ;
        index_mat_3_6(j,4) = uniq_comb_2_vech( uniq_index_3_6(j,3), 0 ) ;
        index_mat_3_6(j,5) = uniq_comb_2_vech( uniq_index_3_6(j,3), 1 ) ;
        index_mat_3_6(j,6) = uniq_comb_2_vech( uniq_index_3_6(j,4), 0 ) ;
        index_mat_3_6(j,7) = uniq_comb_2_vech( uniq_index_3_6(j,4), 1 ) ;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("UniqComb_3rd_4") = join_rows(index_mat_3_4,freq_4) ,
        Rcpp::Named("UniqComb_3rd_6") = join_rows(index_mat_3_6,freq_6) ,
        Rcpp::Named("Uniq3") =uniq_comb_3
    ) ;
}

//MG Laplace


//For GRM, derivative of probbound?
//Check for GRM


//MG Laplace
//Old code
/*
//[[Rcpp::export]]
arma::vec gi(arma::vec theta, arma::vec apars, arma::vec bpars, std::string modeltype, arma::uword mi){
    //  This function compute the probabilities.
    arma::vec out = zeros<vec>(mi);
    if(modeltype == "GPCM"){
        double midouble = mi;
        arma::vec mival = regspace(1.0, midouble);
        arma::vec sumeach = zeros<vec>(mi);
        for(arma::uword cat = 0; cat < mi; cat++){
            arma::uvec ids = find(mival - 1.0 <= cat);
            sumeach(cat) = mival(cat) * sum(apars % theta) + sum(bpars.elem(ids));
        }
        arma::vec expsumeach = exp(sumeach);
        double sumexpsumeach = sum(expsumeach);
        for(arma::uword cat = 0; cat < mi; cat++) out(cat) = expsumeach(cat) / sumexpsumeach;
    }
    if(modeltype == "GRM"){
        arma::vec Pjs = zeros<vec>(mi + 1);
        Pjs(0) = 1.0;
        for(arma::uword cat = 1; cat < mi; cat++){
            double tmp = sum(apars % theta) + bpars(cat);
            tmp = exp(tmp);
            Pjs(cat) = tmp / (1.0 + tmp);
        }
        for(arma::uword cat = 0; cat < mi; cat++){
            out(cat) = Pjs(cat) - Pjs(cat + 1);
        }
    }
    //NRM
	//Need to make apars a matrix or change code below (make vector)
	if(modeltype == "NRM"){
	    arma::vec num = zeros<vec>(mi);
        for(arma::uword cat = 1; cat < mi; cat++) num(cat) = sum(apars.row(cat) % trans(theta)) + bpars(cat);
        arma::vec expnum = exp(num);
        double sumexpnum = sum(expnum);
        for(arma::uword cat = 0; cat < mi; cat++) out(cat) = expnum(cat) / sumexpnum;
	}	
    return(out);
}

//For GRM, derivative of probbound?
//Check for GRM
//[[Rcpp::export]]
arma::mat dgidz(arma::vec apars, std::string modeltype, arma::vec probs, arma::uword mi, arma::uword p){
    arma::mat out = zeros<mat>(p, mi);
    if(modeltype == "GPCM"){
        double midouble = mi;
        arma::vec mival = regspace(1.0, midouble);
        for(arma::uword j = 0; j < p; j++){
            for(arma::uword cat = 0; cat < mi; cat++){
                out(j, cat) = probs(cat) * mival(cat) * apars(j) - probs(cat) * apars(j) * sum(mival % probs.subvec(0, mi-1));
            }
        }
    }
    if(modeltype == "GRM"){
        arma::vec Pis = zeros<vec>(mi + 1);
        Pis(0) = 1.0;
        for(arma::uword cat = 0; cat < mi; cat++) Pis(cat + 1) = Pis(cat) - probs(cat);
        for(arma::uword j = 0; j < p; j++){
            for(arma::uword cat = 0; cat < mi; cat++) out(j, cat) = apars(j) * Pis(cat) * (1.0 - Pis(cat)) - apars(j) * Pis(cat + 1) * (1.0 - Pis(cat + 1));
        }
    }
	if(modeltype == "NRM"){
        for(arma::uword j = 0; j < p; j++){
            for(arma::uword cat = 0; cat < mi; cat++){
                out(j, cat) = probs(cat) * apars(cat, j) - probs(cat) * sum(apars(span(1, mi - 1), span(j)) % probs.subvec(1, mi - 1));
            }
        }
	}
    return(out);
}
*/

//[[Rcpp::export]]
arma::vec gi(arma::vec z, arma::vec apars, arma::vec bpars, std::string modeltype, arma::uword mi, double y){
    //  This function computes the probabilities (for categorical data), the probability mass function (for count data) or the density (for continuous data).
	//  Input is: 
	//  z - a vector with the values of the latent variables
	//  apars - a vector with the discrimination parameters
	//  bpars - a vector with the intercept parameters (for categorical data) or a vector with the intercept parameter and the phi-parameter (for count and continuous data; the phi-parameter is placed after the intercept)
	// modeltype - a string with the measurement model specified ("GPCM", "GRM", "poisson", "normal", "negbin")
	// mi - an integer indicating the number of categories (is 1 for count and continuous data)
	// y - the observed response (for count and continuous data)
    arma::vec out = zeros<vec>(mi);
    if(modeltype == "GPCM"){
        double midouble = mi;
        arma::vec mival = regspace(1.0, midouble);
        arma::vec sumeach = zeros<vec>(mi);
        for(arma::uword cat = 0; cat < mi; cat++){
            arma::uvec ids = find(mival - 1.0 <= cat);
            sumeach(cat) = mival(cat) * sum(apars % z) + sum(bpars.elem(ids));
        }
        arma::vec expsumeach = exp(sumeach);
        double sumexpsumeach = sum(expsumeach);
        for(arma::uword cat = 0; cat < mi; cat++) out(cat) = expsumeach(cat) / sumexpsumeach;
    }
    if(modeltype == "GRM"){
        arma::vec Pjs = zeros<vec>(mi + 1);
        Pjs(0) = 1.0;
        for(arma::uword cat = 1; cat < mi; cat++){
            double tmp = sum(apars % z) + bpars(cat);
            tmp = exp(tmp);
            Pjs(cat) = tmp / (1.0 + tmp);
        }
        for(arma::uword cat = 0; cat < mi; cat++){
            out(cat) = Pjs(cat) - Pjs(cat + 1);
        }
    }
    //NRM
	//Need to make apars a matrix or change code below (make vector)
	if(modeltype == "NRM"){
	    arma::vec num = zeros<vec>(mi);
        for(arma::uword cat = 1; cat < mi; cat++) num(cat) = sum(apars.row(cat) % trans(z)) + bpars(cat);
        arma::vec expnum = exp(num);
        double sumexpnum = sum(expnum);
        for(arma::uword cat = 0; cat < mi; cat++) out(cat) = expnum(cat) / sumexpnum;
	}
	if(modeltype == "negbin"){
		double linpred = bpars(1) + sum(apars % z);
		double phi = bpars(2);
		double topow1 = exp(linpred) / (1.0 / phi + exp(linpred));
		double topow2 = 1.0 / (1.0 + phi * exp(linpred));
		out = (tgamma(y + 1.0 / phi) / (rcpp_factorial(y) * tgamma(1.0 / phi))) * pow(topow1, y) * pow(topow2, 1.0 / phi);
	}
	//Normal
	if(modeltype == "normal"){
		double linpred = bpars(1) + sum(apars % z);
		out = exp((y * linpred - linpred * linpred / 2.0) / bpars(2) - y * y / (2.0 * bpars(2)) - log(2 * M_PI * bpars(2)) / 2.0);
	}
	//Poisson
	if(modeltype == "poisson"){
		double linpred = bpars(1) + sum(apars % z);
		double explinpred = exp(linpred);
		double yfac = rcpp_factorial(y);
		out = pow(explinpred, y) * exp(-explinpred) / yfac;
	}
	//   Output is a vector of probabilities (categorical data) or the probability mass/density (count/continuous data).
    return(out);
}

//[[Rcpp::export]]
arma::mat dgidz(arma::vec z, arma::vec apars, arma::vec bpars, std::string modeltype, arma::vec probs, arma::uword mi, arma::uword p, double y){
	 //  This function computes the first derivatives of a) the probabilities (for categorical data), b) the probability mass function (for count data) or c) the density (for continuous data).
	//  Input is: 
	//  z - a vector with the values of the latent variables
	//  apars - a vector with the discrimination parameters
	//  bpars - a vector with the intercept parameters (for categorical data) or a vector with the intercept parameter and the phi-parameter (for count and continuous data; the phi-parameter is placed after the intercept)
	// modeltype - a string with the measurement model specified ("GPCM", "GRM", "poisson", "normal", "negbin")
	// probs - a vector of probabilities, probability mass function or density at the evaluated latent variables (from function gi())
	// mi - an integer indicating the number of categories (is set to 1 for count and continuous data)
	// p - the number of latent variables related to the observed variable
    arma::mat out = zeros<mat>(p, mi);
    if(modeltype == "GPCM"){
        double midouble = mi;
        arma::vec mival = regspace(1.0, midouble);
        for(arma::uword j = 0; j < p; j++){
            for(arma::uword cat = 0; cat < mi; cat++){
                out(j, cat) = probs(cat) * mival(cat) * apars(j) - probs(cat) * apars(j) * sum(mival % probs.subvec(0, mi-1));
            }
        }
    }
    if(modeltype == "GRM"){
        arma::vec Pis = zeros<vec>(mi + 1);
        Pis(0) = 1.0;
        for(arma::uword cat = 0; cat < mi; cat++) Pis(cat + 1) = Pis(cat) - probs(cat);
        for(arma::uword j = 0; j < p; j++){
            for(arma::uword cat = 0; cat < mi; cat++) out(j, cat) = apars(j) * Pis(cat) * (1.0 - Pis(cat)) - apars(j) * Pis(cat + 1) * (1.0 - Pis(cat + 1));
        }
    }
	if(modeltype == "NRM"){
        for(arma::uword j = 0; j < p; j++){
            for(arma::uword cat = 0; cat < mi; cat++){
                out(j, cat) = probs(cat) * apars(cat, j) - probs(cat) * sum(apars(span(1, mi - 1), span(j)) % probs.subvec(1, mi - 1));
            }
        }
	}
	//if(modeltype == "negbin"){
	//	
	//}
	//Normal
	if(modeltype == "normal"){
		double linpred = bpars(1) + sum(apars % z);
		for(uword j = 0; j < p; j++){
			out(j, 0) = -probs(0) * apars(j) * (linpred - y) / bpars(2);
		}
	}
	//Poisson
	if(modeltype == "poisson"){
		double linpred = bpars(1) + sum(apars % z);
		double explinpred = exp(linpred);
		for(uword j = 0; j < p; j++) {
			out(j, 0) = probs(0) * (y - explinpred) * apars(j);
		}
	}
	//   Output is a matrix of derivatives of the probabilities (categorical data) or the probability mass/density (count/continuous data).
    return(out);
}

// 2023-07-27	Add option for MLE and WLE. MLE is just defining:
//					fgrad = dlh;
//					fhess = d2lh;
// [[Rcpp::export]]
arma::vec optC(arma::vec start, double tol, arma::uword maxit, arma::vec y, arma::uword J, arma::vec cats, arma::uword p, Rcpp::List model, Rcpp::List modelpars, std::vector<std::string> modeltype, std::vector<std::string> link, std::vector<std::string> estimator, arma::vec mu, arma::mat invSigma){
	arma::vec res0 = start;
	arma::uword indexj;
	arma::uword indexk;
	for(arma::uword iter = 0; iter < maxit; iter++){
		Rcpp::checkUserInterrupt();
		arma::vec fgrad;
		arma::mat fhess;
		arma::vec dlh = zeros<vec>(p);
		arma::mat d2lh = zeros<mat>(p, p);
		for(arma::uword i = 0; i < J; i++){
			if(y(i) == 9999) continue;
			Rcpp::List modelparsi = modelpars(i);
			arma::vec apars = modelparsi(0);
			arma::vec bpars = modelparsi(1);
			arma::vec idim = model(i);
			arma::uword pi = idim.n_elem;
			arma::vec thetai = zeros<vec>(pi);
			for(arma::uword j = 0; j < pi; j++){
				indexj = idim(j);
				thetai(j) = res0(indexj);
			}
			//Continuous data: normal
 			if(modeltype[i] == "normal"){
				double ydouble = y(i);
				double linpred = bpars(1) + sum(apars % thetai);
				double phi = bpars(2);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					dlh(indexj) += apars(j) * (linpred - ydouble) / phi;
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						d2lh(indexj, indexk) += apars(j) * apars(k) / phi;
					}
				}
			}
			//Count data: negative binomial
			if(modeltype[i] == "negbin"){
				double ydouble = y(i);
				double linpred = bpars(1) + sum(apars % thetai);
				double explinpred = exp(linpred);
				double phi = bpars(2);
				double topow = (1.0 + phi * explinpred);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					dlh(indexj) += -ydouble * apars(j) + (phi * ydouble + 1.0) * explinpred / topow * apars(j);
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						d2lh(indexj, indexk) += (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * apars(j) * apars(k);
					}
				}
			}
			//Count data: poisson
			if(modeltype[i] == "poisson"){
				double ydouble = y(i);
				double linpred = bpars(1) + sum(apars % thetai);
				double explinpred = exp(linpred);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					dlh(indexj) += (explinpred - ydouble) * apars(j);
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						d2lh(indexj, indexk) += explinpred * apars(j) * apars(k);
					}
				}
			}
			//Ordinal data
			if(modeltype[i] == "GPCM"){
				double ydouble = y(i);
				double midouble = cats(i);
				arma::vec mival = regspace(1.0, midouble);
				arma::vec probi = gi(thetai, apars, bpars, modeltype[i], cats(i), ydouble);
				arma::mat dprobi = dgidz(thetai, apars, bpars, modeltype[i], probi, cats(i), pi, ydouble);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					dlh(indexj) += -apars(j) * (ydouble - sum(probi % mival));
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						d2lh(indexj, indexk) += apars(j) * sum(trans(dprobi.row(k)) % mival);
					}
				}
			}
			if(modeltype[i] == "GRM"){
				if(link[i] == "logit"){
				    arma::uword mi = cats(i);
				    arma::vec Pis = zeros<vec>(mi + 1);
				    Pis(0) = 1.0;
				    for(arma::uword cat = 1; cat < mi; cat++){
				        double tmp = sum(apars % thetai) + bpars(cat);
				        tmp = exp(tmp);
				        Pis(cat) = tmp / (1.0 + tmp); // SJ: I think use plogis() will be better.
				    }
				    for(arma::uword j = 0; j < pi; j++){
				        indexj = idim(j);
				        dlh(indexj) += -apars(j) * (1.0 - Pis(y(i) - 1) - Pis(y(i)));
				        for(arma::uword k = 0; k < pi; k++){
				            indexk = idim(k);
				            d2lh(indexj, indexk) += apars(j) * apars(k) * (Pis(y(i) - 1) * (1.0 - Pis(y(i) - 1)) + Pis(y(i)) * (1.0 - Pis(y(i))));
				        }
				    }
				    
				}
				//else if(link[i] == "probit"){
				    // SJ: Not sure whether this is true...
				    //arma::vec Pi = gj_GRM_probit(y(j), thetaj, apars, bpars, cats(j)) ;//  SJ: Not only Pr().
				    //arma::cube d2Pidt2 = d2gjd2t_GRM_probit(y(n, j), apars, Pi(2), Pi(3), Pi(4), Pi(5), cats(j), pj);
				    //arma::uword index1 = 0;
				    //arma::uword index2;
				    ////This loop does repeated calculations for multidimensional models and can be improved. However, we need to fill out these objects in the end so changing it might not be worth it.
				    //for(auto i : jdim){
				    //    index2 = 0;
				    //    dlh(i) += -apars(index1) * (Pi(4) - Pi(5));
//
				    //    //  SJ: The inner for-loop computes d2h / d2theta (d2hn), the Hessian of h with respect to latent variables
				    //    for(auto j : jdim){
				    //        d2lh(i, j) += -d2Pidt2(index1,index2,0);
				    //        index2 += 1;
				    //    }
				    //    
				    //    index1 += 1;
				    //}
				//}
			}
			if(modeltype[i] == "NRM"){
				arma::uword mi = cats(i);
				double ydouble = y(i);
				arma::vec probi = gi(thetai, apars, bpars, modeltype[i], mi, ydouble);
				arma::mat dprobi = dgidz(thetai, apars, bpars, modeltype[i], probi, cats(i), pi, ydouble);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					dlh(indexj) += -(ydouble - sum(apars(span(1, mi - 1), span(j)) % probi.subvec(1, mi - 1)));
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						d2lh(indexj, indexk) += sum(apars(span(1, mi - 1), span(j)) % trans(dprobi.row(k)));
					}
				}
			}
		}
		if(estimator[0] == "MLE"){
			fgrad = dlh;
			fhess = d2lh;
		} else if(estimator[0] == "MAP"){
			fgrad = trans(invSigma) * (res0 - mu) + dlh;
			fhess = d2lh + invSigma;				
		} 
		//else if(estimator == "WLE"){
		//	fgrad = ;
		//	fhess = ;				
		//}
		//fgrad = trans(invSigma) * (res0 - mu) + dlh;
		//fhess = d2lh + invSigma;
		arma::vec res1 = res0 - solve(fhess, fgrad);
		if(all(abs((res0 - res1) / res0) < tol)) return(res1);
		res0 = res1;
	}
	arma::vec toret = zeros<vec>(p);
	toret(0) = 99.0;
	return(toret);
}


// [[Rcpp::export]]
arma::mat infoC(arma::vec start, double tol, arma::uword maxit, arma::vec y, arma::uword J, arma::vec cats, arma::uword p, Rcpp::List model, Rcpp::List modelpars, std::vector<std::string> modeltype, std::vector<std::string> link, std::vector<std::string> estimator, std::vector<std::string> information, arma::vec mu, arma::mat invSigma){
	arma::vec res0 = start;
	arma::uword indexj;
	arma::uword indexk;
	arma::mat fhess;
	//for(arma::uword iter = 0; iter < maxit; iter++){
		Rcpp::checkUserInterrupt();
		//arma::vec fgrad;
		//arma::vec dlh = zeros<vec>(p);
		arma::mat d2lh = zeros<mat>(p, p);
		for(arma::uword i = 0; i < J; i++){
			if(y(i) == 9999) continue;
			Rcpp::List modelparsi = modelpars(i);
			arma::vec apars = modelparsi(0);
			arma::vec bpars = modelparsi(1);
			arma::vec idim = model(i);
			arma::uword pi = idim.n_elem;
			arma::vec thetai = zeros<vec>(pi);
			for(arma::uword j = 0; j < pi; j++){
				indexj = idim(j);
				thetai(j) = res0(indexj);
			}
			//Continuous data: normal
 			if(modeltype[i] == "normal"){
				//double ydouble = y(i);
				//double linpred = bpars(1) + sum(apars % thetai);
				double phi = bpars(2);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					//dlh(indexj) += apars(j) * (linpred - ydouble) / phi;
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						if(information[0] == "observed"){
							d2lh(indexj, indexk) += apars(j) * apars(k) / phi;
						} else if(information[0] == "expected"){
							d2lh(indexj, indexk) += apars(j) * apars(k) / phi;
						}
					}
				}
			}
			//Count data: negative binomial
			if(modeltype[i] == "negbin"){
				double ydouble = y(i);
				double linpred = bpars(1) + sum(apars % thetai);
				double explinpred = exp(linpred);
				double phi = bpars(2);
				double topow = (1.0 + phi * explinpred);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					//dlh(indexj) += -ydouble * apars(j) + (phi * ydouble + 1.0) * explinpred / topow * apars(j);
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						if(information[0] == "observed"){
							d2lh(indexj, indexk) += (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * apars(j) * apars(k);
						} else if(information[0] == "expected"){
							//Update here
							d2lh(indexj, indexk) += explinpred / topow * apars(j) * apars(k);
						}
					}
				}
			}
			//Count data: poisson
			if(modeltype[i] == "poisson"){
				//double ydouble = y(i);
				double linpred = bpars(1) + sum(apars % thetai);
				double explinpred = exp(linpred);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					//dlh(indexj) += (explinpred - ydouble) * apars(j);
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						if(information[0] == "observed"){
							d2lh(indexj, indexk) += explinpred * apars(j) * apars(k);
						} else if(information[0] == "expected"){
							d2lh(indexj, indexk) += explinpred * apars(j) * apars(k);
						}
					}
				}
			}
			//Ordinal data
			if(modeltype[i] == "GPCM"){
				double ydouble = y(i);
				double midouble = cats(i);
				arma::vec mival = regspace(1.0, midouble);
				arma::vec probi = gi(thetai, apars, bpars, modeltype[i], cats(i), ydouble);
				arma::mat dprobi = dgidz(thetai, apars, bpars, modeltype[i], probi, cats(i), pi, ydouble);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					//dlh(indexj) += -apars(j) * (ydouble - sum(probi % mival));
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						if(information[0] == "observed"){
							d2lh(indexj, indexk) += apars(j) * sum(trans(dprobi.row(k)) % mival);
						} else if(information[0] == "expected"){
							d2lh(indexj, indexk) += apars(j) * sum(trans(dprobi.row(k)) % mival);
						}
					}
				}
			}
			if(modeltype[i] == "GRM"){
				if(link[i] == "logit"){
				    arma::uword mi = cats(i);
				    arma::vec Pis = zeros<vec>(mi + 1);
					arma::vec Pi = zeros<vec>(mi);
					double tempinfo;
				    Pis(0) = 1.0;
				    for(arma::uword cat = 1; cat < mi; cat++){
				        double tmp = sum(apars % thetai) + bpars(cat);
				        tmp = exp(tmp);
				        Pis(cat) = tmp / (1.0 + tmp); // SJ: I think use plogis() will be better.
				    }
					if(information[0] == "expected"){
						for(arma::uword cat = 0; cat < mi; cat++){
							Pi(cat) = Pis(cat) - Pis(cat + 1);
						}
					}
				    for(arma::uword j = 0; j < pi; j++){
				        indexj = idim(j);
				        //dlh(indexj) += -apars(j) * (1.0 - Pis(y(i) - 1) - Pis(y(i)));
				        for(arma::uword k = 0; k < pi; k++){
				            indexk = idim(k);
							if(information[0] == "observed"){
								d2lh(indexj, indexk) += apars(j) * apars(k) * (Pis(y(i) - 1) * (1.0 - Pis(y(i) - 1)) + Pis(y(i)) * (1.0 - Pis(y(i))));
							} else if(information[0] == "expected"){
								tempinfo = 0.0;
								for(arma::uword cat = 0; cat < mi; cat++){
									tempinfo += Pi(cat) * apars(j) * apars(k) * (Pis(cat) * (1.0 - Pis(cat)) + Pis(cat + 1) * (1.0 - Pis(cat + 1)));
								}								
								d2lh(indexj, indexk) += tempinfo;
							}
				        }
				    }
				    
				}
				//else if(link[i] == "probit"){
				    // SJ: Not sure whether this is true...
				    //arma::vec Pi = gj_GRM_probit(y(j), thetaj, apars, bpars, cats(j)) ;//  SJ: Not only Pr().
				    //arma::cube d2Pidt2 = d2gjd2t_GRM_probit(y(n, j), apars, Pi(2), Pi(3), Pi(4), Pi(5), cats(j), pj);
				    //arma::uword index1 = 0;
				    //arma::uword index2;
				    ////This loop does repeated calculations for multidimensional models and can be improved. However, we need to fill out these objects in the end so changing it might not be worth it.
				    //for(auto i : jdim){
				    //    index2 = 0;
				    //    dlh(i) += -apars(index1) * (Pi(4) - Pi(5));
//
				    //    //  SJ: The inner for-loop computes d2h / d2theta (d2hn), the Hessian of h with respect to latent variables
				    //    for(auto j : jdim){
				    //        d2lh(i, j) += -d2Pidt2(index1,index2,0);
				    //        index2 += 1;
				    //    }
				    //    
				    //    index1 += 1;
				    //}
				//}
			}
			if(modeltype[i] == "NRM"){
				arma::uword mi = cats(i);
				double ydouble = y(i);
				arma::vec probi = gi(thetai, apars, bpars, modeltype[i], mi, ydouble);
				arma::mat dprobi = dgidz(thetai, apars, bpars, modeltype[i], probi, cats(i), pi, ydouble);
				for(arma::uword j = 0; j < pi; j++){
					indexj = idim(j);
					//dlh(indexj) += -(ydouble - sum(apars(span(1, mi - 1), span(j)) % probi.subvec(1, mi - 1)));
					for(arma::uword k = 0; k < pi; k++){
						indexk = idim(k);
						if(information[0] == "observed"){
							d2lh(indexj, indexk) += sum(apars(span(1, mi - 1), span(j)) % trans(dprobi.row(k)));
						} else if(information[0] == "expected"){
							d2lh(indexj, indexk) += sum(apars(span(1, mi - 1), span(j)) % trans(dprobi.row(k)));
						}
					}
				}
			}
		}
		if(estimator[0] == "MLE"){
			//fgrad = dlh;
			fhess = d2lh;
		} else if(estimator[0] == "MAP"){
			//fgrad = trans(invSigma) * (res0 - mu) + dlh;
			fhess = d2lh + invSigma;				
		} 
		//else if(estimator == "WLE"){
		//	fgrad = ;
		//	fhess = ;				
		//}
		//fgrad = trans(invSigma) * (res0 - mu) + dlh;
		//fhess = d2lh + invSigma;
		//arma::vec res1 = res0 - solve(fhess, fgrad);
		//if(all(abs((res0 - res1) / res0) < tol)) return(res1);
		//res0 = res1;
	//}
	//arma::vec toret = zeros<vec>(p);
	//toret(0) = 99.0;
	//return(toret);
	return(fhess);
}


// [[Rcpp::export]]
double hoptC(arma::vec theta, arma::vec y, arma::uword J, arma::vec cats, arma::uword p, Rcpp::List model, Rcpp::List modelpars, std::vector<std::string> modeltype, std::vector<std::string> link, arma::vec mu, arma::mat Sigma){
    double lh = 0.0;
    arma::mat tmpsig = 2.0 * M_PI * Sigma;
    arma::mat invsigma = inv(Sigma);
    double logsqrtdet2pisigma = det(tmpsig);
    logsqrtdet2pisigma = sqrt(logsqrtdet2pisigma);
    logsqrtdet2pisigma = log(logsqrtdet2pisigma);
    for(arma::uword i = 0; i < J; i++){
        arma::uword yi = 1;
		if(modeltype[i] == "GPCM"){
			yi = y(i);
		} else if(modeltype[i] == "GRM"){
			yi = y(i);
		} else if(modeltype[i] == "negbin"){
			yi = 1;
		} else if(modeltype[i] == "normal"){
			yi = 1;
		} else if(modeltype[i] == "poisson"){
			yi = 1;
		}
		double ydouble = y(i);
        double probb;
        if(y(i) == 9999){
            probb = 1.0;
            lh += log(probb);
            continue;
        }
        Rcpp::List modelparsi = modelpars(i);
        arma::vec apars = modelparsi(0);
        arma::vec bpars = modelparsi(1);
        arma::vec idim = model(i);
        arma::uword pi = idim.n_elem;
        arma::vec thetai = zeros<vec>(pi);
        for(arma::uword j = 0; j < pi; j++){
            arma::uword indexj = idim(j);
            thetai(j) = theta(indexj);
        }
        arma::uword cati = cats(i);
        arma::vec probj = gi(thetai, apars, bpars, modeltype[i], cati, ydouble);
        probb = probj(yi - 1);
        lh += log(probb);
    }
    arma::mat tmpp = -0.5 * trans(theta - mu) * invsigma * (theta - mu);
    lh = tmpp(0,0) - logsqrtdet2pisigma + lh;	
    return(-lh);
}

// [[Rcpp::export]]
arma::vec dhoptC(arma::vec theta, arma::vec y, arma::uword J, arma::vec cats, arma::uword p, Rcpp::List model, Rcpp::List modelpars, std::vector<std::string> modeltype, std::vector<std::string> link, arma::vec mu, arma::mat Sigma){
    arma::vec dlh = zeros<vec>(p);
    for(arma::uword i = 0; i < J; i++){
        arma::uword yi = 1;
		if(modeltype[i] == "GPCM"){
			yi = y(i);
		} else if(modeltype[i] == "GRM"){
			yi = y(i);
		} else if(modeltype[i] == "negbin"){
			yi = 1;
		} else if(modeltype[i] == "normal"){
			yi = 1;
		} else if(modeltype[i] == "poisson"){
			yi = 1;
		}
		double ydouble = y(i);
        if(y(i) == 9999) continue;
        Rcpp::List modelparsi = modelpars(i);
        arma::vec apars = modelparsi(0);
        arma::vec bpars = modelparsi(1);
        arma::vec idim = model(i);
        arma::uword pi = idim.n_elem;
        arma::vec thetai = zeros<vec>(pi);
        for(arma::uword j = 0; j < pi; j++){
            arma::uword indexj = idim(j);
            thetai(j) = theta(indexj);
        }
		if(modeltype[i] == "negbin"){
			double ydouble = y(i);
			double linpred = bpars(1) + sum(apars % thetai);
			double explinpred = exp(linpred);
			double phi = bpars(2);
			double topow = (1.0 + phi * explinpred);
			for(arma::uword j = 0; j < pi; j++){
				arma::uword indexj = idim(j);
				dlh(indexj) += -ydouble * apars(j) + (phi * ydouble + 1.0) * explinpred / topow * apars(j);
			}
		} else{
			arma::vec probi = gi(thetai, apars, bpars, modeltype[i], cats(i), ydouble);
			arma::mat dprobi = dgidz(thetai, apars, bpars, modeltype[i], probi, cats(i), pi, ydouble);
			for(arma::uword j = 0; j < pi; j++){
				arma::uword indexj = idim(j);
				double probb = probi(yi - 1);
				dlh(indexj) += -(1.0 / probb) * dprobi(j, yi - 1);
			}
		}
    }
    arma::mat invsigma = inv(Sigma);
    dlh = trans(invsigma) * (theta - mu) + dlh;
    return(dlh);
}

// [[Rcpp::export]]
double mleoptC(arma::vec theta, arma::vec y, arma::uword J, arma::vec cats, arma::uword p, Rcpp::List model, Rcpp::List modelpars, std::vector<std::string> modeltype, std::vector<std::string> link){
    double lh = 0.0;
    for(arma::uword i = 0; i < J; i++){
        arma::uword yi = 1;
		if(modeltype[i] == "GPCM"){
			yi = y(i);
		} else if(modeltype[i] == "GRM"){
			yi = y(i);
		} else if(modeltype[i] == "negbin"){
			yi = 1;
		} else if(modeltype[i] == "normal"){
			yi = 1;
		} else if(modeltype[i] == "poisson"){
			yi = 1;
		}
		double ydouble = y(i);
        double probb;
        if(y(i) == 9999){
            probb = 1.0;
            lh += log(probb);
            continue;
        }
        Rcpp::List modelparsi = modelpars(i);
        arma::vec apars = modelparsi(0);
        arma::vec bpars = modelparsi(1);
        arma::vec idim = model(i);
        arma::uword pi = idim.n_elem;
        arma::vec thetai = zeros<vec>(pi);
        for(arma::uword j = 0; j < pi; j++){
            arma::uword indexj = idim(j);
            thetai(j) = theta(indexj);
        }
        arma::uword cati = cats(i);
        arma::vec probj = gi(thetai, apars, bpars, modeltype[i], cati, ydouble);
        probb = probj(yi - 1);
        lh += log(probb);
    }
    return(-lh);
}

// [[Rcpp::export]]
arma::vec dmleoptC(arma::vec theta, arma::vec y, arma::uword J, arma::vec cats, arma::uword p, Rcpp::List model, Rcpp::List modelpars, std::vector<std::string> modeltype, std::vector<std::string> link){
    arma::vec dlh = zeros<vec>(p);
    for(arma::uword i = 0; i < J; i++){
        arma::uword yi = 1;
		if(modeltype[i] == "GPCM"){
			yi = y(i);
		} else if(modeltype[i] == "GRM"){
			yi = y(i);
		} else if(modeltype[i] == "negbin"){
			yi = 1;
		} else if(modeltype[i] == "normal"){
			yi = 1;
		} else if(modeltype[i] == "poisson"){
			yi = 1;
		}
		double ydouble = y(i);
        if(y(i) == 9999) continue;
        Rcpp::List modelparsi = modelpars(i);
        arma::vec apars = modelparsi(0);
        arma::vec bpars = modelparsi(1);
        arma::vec idim = model(i);
        arma::uword pi = idim.n_elem;
        arma::vec thetai = zeros<vec>(pi);
        for(arma::uword j = 0; j < pi; j++){
            arma::uword indexj = idim(j);
            thetai(j) = theta(indexj);
        }
		if(modeltype[i] == "negbin"){
			double ydouble = y(i);
			double linpred = bpars(1) + sum(apars % thetai);
			double explinpred = exp(linpred);
			double phi = bpars(2);
			double topow = (1.0 + phi * explinpred);
			for(arma::uword j = 0; j < pi; j++){
				arma::uword indexj = idim(j);
				dlh(indexj) += -ydouble * apars(j) + (phi * ydouble + 1.0) * explinpred / topow * apars(j);
			}
		} else{
			arma::vec probi = gi(thetai, apars, bpars, modeltype[i], cats(i), ydouble);
			arma::mat dprobi = dgidz(thetai, apars, bpars, modeltype[i], probi, cats(i), pi, ydouble);
			for(arma::uword j = 0; j < pi; j++){
				arma::uword indexj = idim(j);
				double probb = probi(yi - 1);
				dlh(indexj) += -(1.0 / probb) * dprobi(j, yi - 1);
			}
		}
    }
    return(dlh);
}


// [[Rcpp::export]]
arma::mat dhdb(arma::mat theta, arma::mat sigma, arma::mat betamat, arma::mat x, arma::uword ndim, arma::uword nind){
	arma::uword nx = x.n_cols - 1;
	arma::uword ntotpar = ndim * nx;
	arma::mat out = zeros<mat>(nind, ntotpar);
	arma::vec meanvec = zeros<mat>(nx);
	for(arma::uword k = 0; k < nx; k++) meanvec(k) = mean(x.col(k + 1));
	arma::mat invsigma = inv(sigma);
	for(arma::uword i = 0; i < nind; i++){
		arma::vec tempvec = trans(x.row(i));
		arma::vec constvec = zeros<vec>(ndim);
		arma::vec tempvec1 = trans(theta.row(i)) - betamat * tempvec;
		constvec = invsigma * tempvec1;
		arma::uword tmp = 0;
		for(arma::uword j = 0; j < ndim; j++){
			for(arma::uword k = 0; k < nx; k++){
				out(i, tmp) = constvec(j) * meanvec(k) - constvec(j) * tempvec(k + 1);
				tmp = tmp + 1;
			}
		}
	}
	return(out);
}

// [[Rcpp::export]]
arma::cube d2hdtdb(arma::mat theta, arma::mat sigma, arma::mat betamat, arma::mat x, arma::uword ndim, arma::uword nind){
	arma::uword nx = x.n_cols - 1;
	arma::uword ntotpar = ndim * nx;
	arma::cube out = zeros<cube>(nind, ndim, ntotpar);
	arma::vec meanvec = zeros<vec>(nx);
	for(arma::uword k = 0; k < nx; k++) meanvec(k) = mean(x.col(k + 1));
	arma::mat invsigma = inv(sigma);
	for(arma::uword i = 0; i < nind; i++){
		arma::uword tmp = 0;
		for(arma::uword j = 0; j < ndim; j++){
			for(arma::uword k = 0; k < nx; k++){
				for(arma::uword jj = 0; jj < ndim; jj++) out(i, jj, tmp) = invsigma(jj, j) * (-meanvec(k) + x(i, k + 1));
				tmp = tmp + 1;
			}
		}
	}
	return(out);
}

//Old code
/*
//[[Rcpp::export]]
arma::cube d2gidz2(arma::vec apars, std::string modeltype, arma::vec probs, arma::mat dprobs, arma::uword mi, arma::uword p){
    arma::cube out = zeros<cube>(p, p, mi);
    double midouble = mi;
    arma::vec mival = regspace(1.0, midouble);
    for(arma::uword j = 0; j < p; j++){
        for(arma::uword k = 0; k < p; k++){
            for(arma::uword cat = 0; cat < mi; cat++){
                out(j, k, cat) = dprobs(k, cat) * mival(cat) * apars(j) - dprobs(k, cat) * apars(j) * sum(mival % probs) - probs(cat) * apars(j) * sum(trans(dprobs.row(k)) % mival);
            }
        }
    }
    return(out);
}

arma::mat dgidu(arma::vec theta, std::string modeltype, arma::vec probs, arma::uword mi, arma::uword p){
    arma::mat out;
    arma::vec sumeach = zeros<vec>(mi);
    double midouble = mi;
    arma::vec mival = regspace(1.0, midouble);
    if(modeltype == "GPCM"){
        arma::uword npar = p + mi - 1;
        out = zeros<mat>(npar, mi);
        for(arma::uword cat = 0; cat < mi; cat++){
            //apars
            for(arma::uword u = 0; u < p; u++){
                out(u, cat) = probs(cat) * mival(cat) * theta(u) - probs(cat) * theta(u) * sum(mival % probs);
            }
            //bpars
            arma::uword uu = 1;
            for(arma::uword u = p; u < npar; u++){
                arma::uvec ids = find(mival - 1.0 >= uu);
                if(cat >= uu) out(u, cat) = probs(cat) - probs(cat) * sum(probs.elem(ids));
                if(cat < uu) out(u, cat) = - probs(cat) * sum(probs.elem(ids));
                uu++;
            }
        }
    }
    return(out);	
}

arma::cube d2gidzdu(arma::vec apars, std::string modeltype, arma::vec probs, arma::mat dprobs, arma::mat dprobsdu, arma::uword mi, arma::uword p){
    arma::cube out;
    arma::vec sumeach = zeros<vec>(mi);
    double midouble = mi;
    arma::vec mival = regspace(1.0, midouble);
    if(modeltype == "GPCM"){
        arma::uword npar = p + mi - 1;
        out = zeros<cube>(p, npar, mi);
        for(arma::uword j = 0; j < p; j++){
            for(arma::uword u = 0; u < npar; u++){
                double myind = j == u;
                for(arma::uword cat = 0; cat < mi; cat++){
                    out(j, u, cat) = myind * dprobs(j, cat) / apars(j) + dprobsdu(u, cat) * mival(cat) * apars(j) - dprobsdu(u, cat) * apars(j) * sum(probs % mival) - probs(cat) * apars(j) * sum(trans(dprobsdu.row(u)) % mival);
                }
            }
        }
    }
    return(out);
}
*/

//[[Rcpp::export]]
arma::cube d2gidz2(arma::vec z, arma::vec apars, arma::vec bpars, std::string modeltype, arma::vec probs, arma::mat dprobs, arma::uword mi, arma::uword p, double y){
	//  This function computes the second derivatives of a) the probabilities (for categorical data), b) the probability mass function (for count data) or c) the density (for continuous data).
	//  Input is: 
	//  z - a vector with the values of the latent variables
	//  apars - a vector with the discrimination parameters
	//  bpars - a vector with the intercept parameters (for categorical data) or a vector with the intercept parameter and the phi-parameter (for count and continuous data; the phi-parameter is placed after the intercept)
	// modeltype - a string with the measurement model specified ("GPCM", "GRM", "poisson", "normal", "negbin")
	// probs - a vector of probabilities, probability mass function or density at the evaluated latent variables (from function gi())
	// dprobs - a matrix of first derivatives of probabilities, probability mass function or density at the evaluated latent variables (from function dgidz())
	// mi - an integer indicating the number of categories (is 1 for count and continuous data)
	// p - the number of latent variables
    arma::cube out = zeros<cube>(p, p, mi);
	if(modeltype == "GPCM"){
		double midouble = mi;
		arma::vec mival = regspace(1.0, midouble);
		for(arma::uword j = 0; j < p; j++){
			for(arma::uword k = 0; k < p; k++){
				for(arma::uword cat = 0; cat < mi; cat++){
					out(j, k, cat) = dprobs(k, cat) * mival(cat) * apars(j) - dprobs(k, cat) * apars(j) * sum(mival % probs) - probs(cat) * apars(j) * sum(trans(dprobs.row(k)) % mival);
				}
			}
		}
	}
	//if(modeltype == "negbin"){
	//	
	//}
	//Normal
	if(modeltype == "normal"){
		double linpred = bpars(1) + sum(apars % z);
		for(arma::uword j = 0; j < p; j++) {
			for(arma::uword k = 0; k < p; k++) {
				out(j,k,0) = -apars(j) / bpars(2) * (probs(0) * apars(k) + dprobs(k,0) * (linpred - y));
			}
		}
	}
	//Poisson
	if(modeltype == "poisson"){
		double linpred = bpars(1) + sum(apars % z);
		double explinpred = exp(linpred);
		for(arma::uword j = 0; j < p; j++) {
			for(arma::uword k = 0; k < p; k++) {
				out(j,k,0) = -apars(j) * apars(k) * probs(0) * explinpred + dprobs(k,0) * apars(j) * (y - explinpred);
			}
		}
	}
	//   Output is a cube of second derivatives of the probabilities (categorical data) or the probability mass/density (count/continuous data).
    return(out);
}

//[[Rcpp::export]]
arma::mat dgidu(arma::vec z, arma::vec apars, arma::vec bpars, std::string modeltype, arma::vec probs, arma::uword mi, arma::uword p, double y){
	//  This function computes the first derivatives of a) the probabilities (for categorical data), b) the probability mass function (for count data) or c) the density (for continuous data), with respect to the unknown parameters.
	//  Input is: 
	//  z - a vector with the values of the latent variables
	//  apars - a vector with the discrimination parameters
	//  bpars - a vector with the intercept parameters (for categorical data) or a vector with the intercept parameter and the phi-parameter (for count and continuous data; the phi-parameter is placed after the intercept)
	//  modeltype - a string with the measurement model specified ("GPCM", "GRM", "poisson", "normal", "negbin")
	//  probs - a vector of probabilities, probability mass function or density at the evaluated latent variables (from function gi())
	//  mi - an integer indicating the number of categories (is 1 for count and continuous data)
	//  p - the number of latent variables
    arma::mat out;
    arma::vec sumeach = zeros<vec>(mi);
    double midouble = mi;
    arma::vec mival = regspace(1.0, midouble);
    if(modeltype == "GPCM"){
        arma::uword npar = p + mi - 1;
        out = zeros<mat>(npar, mi);
        for(arma::uword cat = 0; cat < mi; cat++){
            //apars
            for(arma::uword u = 0; u < p; u++){
                out(u, cat) = probs(cat) * mival(cat) * z(u) - probs(cat) * z(u) * sum(mival % probs);
            }
            //bpars
            arma::uword uu = 1;
            for(arma::uword u = p; u < npar; u++){
                arma::uvec ids = find(mival - 1.0 >= uu);
                if(cat >= uu) out(u, cat) = probs(cat) - probs(cat) * sum(probs.elem(ids));
                if(cat < uu) out(u, cat) = - probs(cat) * sum(probs.elem(ids));
                uu++;
            }
        }
    }
	//if(modeltype == "negbin"){
	//	
	//}
	//Normal
	if(modeltype == "normal"){
		arma::uword npar = p + 2;
		out = zeros<mat>(npar,mi);
		double linpred = bpars(1) + sum(apars % z);
		//apars
		for(arma::uword u = 0; u < p; u++){
			out(u,0) = (y - linpred) * probs(0) * z(u) / bpars(2) ;
		}
		//bpars
		out(p,0) = (y - linpred) * probs(0) / bpars(2);
		//phipars
		out(p+1,0) = probs(0) / (2.0 * bpars(2)) * (pow((y - linpred), 2.0) / bpars(2) - 1 );
	}
	//Poisson
	if(modeltype == "poisson"){
		arma::uword npar = p + 1;
		out = zeros<mat>(npar,mi);
		double linpred = bpars(1) + sum(apars % z);
		double explinpred = exp(linpred);
		//apars
		for(arma::uword u = 0; u < p; u++){
			out(u,0) = probs(0) * z(u) * (y - explinpred);
		}
		//bpars
		out(p,0) = probs(0) * (y - explinpred);
	}
	//   Output is a matrix of derivatives of the probabilities (categorical data) or the probability mass/density (count/continuous data).
    return(out);	
}

//[[Rcpp::export]]
arma::cube d2gidzdu(arma::vec z, arma::vec apars, arma::vec bpars, std::string modeltype, arma::vec probs, arma::mat dprobs, arma::mat dprobsdu, arma::uword mi, arma::uword p, double y){
	//  This function computes the first derivatives of first derivatives of a) the probabilities (for categorical data), b) the probability mass function (for count data) or c) the density (for continuous data), with respect to the unknown parameters.
	//  Input is: 
	//  z - a vector with the values of the latent variables
	//  apars - a vector with the discrimination parameters
	//  bpars - a vector with the intercept parameters (for categorical data) or a vector with the intercept parameter and the phi-parameter (for count and continuous data; the phi-parameter is placed after the intercept)
	// modeltype - a string with the measurement model specified ("GPCM", "GRM", "poisson", "normal", "negbin")
	// probs - a vector of probabilities, probability mass function or density at the evaluated latent variables (from function gi())
	// dprobs - a matrix of first derivatives of probabilities, probability mass function or density at the evaluated latent variables (from function dgidz())
	// dprobsdu - a matrix of first derivatives of probabilities, probability mass function or density at the evaluated latent variables (from function dgidu())
	// mi - an integer indicating the number of categories (is 1 for count and continuous data)
	// p - the number of latent variables
    arma::cube out;
    arma::vec sumeach = zeros<vec>(mi);
    double midouble = mi;
    arma::vec mival = regspace(1.0, midouble);
    if(modeltype == "GPCM"){
        arma::uword npar = p + mi - 1;
        out = zeros<cube>(p, npar, mi);
        for(arma::uword j = 0; j < p; j++){
            for(arma::uword u = 0; u < npar; u++){
                double myind = j == u;
                for(arma::uword cat = 0; cat < mi; cat++){
                    out(j, u, cat) = myind * dprobs(j, cat) / apars(j) + dprobsdu(u, cat) * mival(cat) * apars(j) - dprobsdu(u, cat) * apars(j) * sum(probs % mival) - probs(cat) * apars(j) * sum(trans(dprobsdu.row(u)) % mival);
                }
            }
        }
    }
	//if(modeltype == "negbin"){
	//	
	//}
	//Normal
	if(modeltype == "normal"){
		arma::uword npar = p + 2;
		out = zeros<cube>(p, npar, mi);
		double linpred = bpars(1) + sum(apars % z);
		for(arma::uword j = 0; j < p; j++){
			//apars
			for(arma::uword u = 0; u < p; u++){
				double myind = j==u;
				out(j, u, 0) = apars(j) / bpars(2) * ( dprobsdu(u) * (y - linpred) - probs(0) * z(u)) + myind * (y - linpred) * probs(0) / bpars(2);
			}
			//bpars
			out(j, p, 0) = apars(j) / bpars(2) * ( dprobsdu(p) * (y - linpred) - probs(0));
			//phipars
			out(j, p + 1, 0) = apars(j) / bpars(2) * (y - linpred) * (dprobsdu(p + 1) - probs(0) / bpars(2));
		}
	}
	//Poisson
	if(modeltype == "poisson"){
		arma::uword npar = p + 1;
		out = zeros<cube>(p, npar, mi);
		double linpred = bpars(1) + sum(apars % z);
		double explinpred = exp(linpred);
		for(arma::uword j = 0; j < p; j++){
			//apars
			for(arma::uword u = 0; u < p; u++){
				double myind = j==u;
				out(j, u, 0) = dprobsdu(u) * (y - explinpred) * apars(j) - explinpred * probs(0) * z(u) * apars(j) + myind * (y - explinpred) * probs(0);
			}
			//bpars
			out(j, p, 0) = apars(j) * (dprobsdu(p) * (y - explinpred) - probs(0) * explinpred);
		}
	}
	//   Output is a cube of derivatives of the derivatives of the probabilities (categorical data) or the probability mass/density (count/continuous data).
    return(out);
}

// [[Rcpp::export]]
arma::vec gj_GRM_probit(arma::uword ynj, arma::vec theta, arma::vec apars, arma::vec bpars, arma::uword mj){
    //  This function compute the density, probability, and others needed in NRM.
    //  We return ( Pr(1|theta), Pr(2|theta), linear predictor (1), linear predictor (2), ratio (1), ratio (2) ).
    //  The ratio dnorm()/pnorm() is not always needed, perhaps an if-else to determine whether we need it.
    arma::vec out = zeros<vec>(6);
    double sumat = sum(apars % theta);
    if(ynj == 1){
        out(3) = -sumat - bpars(ynj);
        out(0) = 1.0;
        out(1) = R::pnorm( out(3), 0.0, 1.0, 0L, 0L ); //  lower=false.
        out(5) = R::dnorm( out(3), 0.0, 1.0, 0L ) / (out(0) - out(1)); 
    }
    else if(ynj == mj){
        out(2) = -sumat - bpars(ynj - 1);
        out(0) = R::pnorm( out(2), 0.0, 1.0, 0L, 0L );
        out(4) = R::dnorm( out(2), 0.0, 1.0, 0L ) / (out(0) - out(1)); 
    }
    else{
        out(2) = -sumat - bpars(ynj - 1);
        out(3) = -sumat - bpars(ynj);
        out(0) = R::pnorm( out(2), 0.0, 1.0, 0L, 0L ); 
        out(1) = R::pnorm( out(3), 0.0, 1.0, 0L, 0L );
        out(4) = R::dnorm( out(2), 0.0, 1.0, 0L ) / (out(0) - out(1)); 
        out(5) = R::dnorm( out(3), 0.0, 1.0, 0L ) / (out(0) - out(1)); 
    }
    return(out);
}

// [[Rcpp::export]]
arma::mat dgjdu_GRM_probit(arma::uword ynj, arma::vec theta, double dpnorm1, double dpnorm2, arma::uword mj, arma::uword pj, arma::uword nparj){
    //  This function computes dlogPr / dpar.
    arma::mat dPidu = zeros<mat>(nparj, 2);
    //We need the entries that correspond to ynj and ynj + 1.
    if(ynj == 1){
        //a-par
        for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 1) = dpnorm2 * theta(temppar);
        //b-par
        dPidu(pj, 1) = dpnorm2;
    }
    else if(ynj == mj){
        //a-par
        for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 0) = dpnorm1 * theta(temppar);
        //b-par
        dPidu(ynj - 2 + pj, 0) = dpnorm1;
    }
    else{
        //a-par
        for(arma::uword temppar = 0; temppar < pj; temppar++){
            dPidu(temppar, 0) = dpnorm1 * theta(temppar);
            dPidu(temppar, 1) = dpnorm2 * theta(temppar);
        }
        //b-par
        dPidu(ynj - 2 + pj, 0) = dpnorm1;
        dPidu(ynj - 1 + pj, 1) = dpnorm2;
    }
    return(dPidu);
}

// [[Rcpp::export]]
arma::mat dgjdt_GRM_probit(arma::vec apars, double dpnorm1, double dpnorm2, arma::uword pj){
    //  SJ: Compute dlogPi / dtheta
    arma::mat dlogPidt = zeros<mat>(pj, 2);
    for(arma::uword uu = 0; uu < pj; uu++){
        dlogPidt(uu, 0) = dpnorm1 * apars(uu) ;
        dlogPidt(uu, 1) = dpnorm2 * apars(uu) ;
    }
    return(dlogPidt);
}

// [[Rcpp::export]]
arma::cube d2gjd2t_GRM_probit(arma::uword ynj, arma::vec apars, double linpred1, double linpred2, double dpnorm1, double dpnorm2, arma::uword mj, arma::uword pj){
    //  SJ: Compute d2logPi / d2theta
    arma::cube d2logPid2t = zeros<cube>(pj, pj, 1);
    arma::mat aparmat = conv_to< arma::mat >::from(apars); //  This is a column matrix.
    //Rcout << "nrow=" << aparmat.n_rows << ", ncol=" << aparmat.n_cols << std::endl;
    if(ynj == 1){
        d2logPid2t.slice(0) = -1.0*(linpred2 * dpnorm2 + pow(dpnorm2,2.0)) * aparmat * aparmat.t();
    }
    else if(ynj == mj){
        d2logPid2t.slice(0) = (linpred1 * dpnorm1 - pow(dpnorm1,2.0)) * aparmat * aparmat.t();
    }
    else {
        //  SJ: This actually covers both ynj=1 and ynj=mj.
        d2logPid2t.slice(0) = (linpred1 * dpnorm1 - linpred2 * dpnorm2 - pow(dpnorm1 - dpnorm2,2.0)) * aparmat * aparmat.t();
    }
    return(d2logPid2t);
}

arma::cube d2gjdtdu_GRM_probit(arma::uword ynj, arma::vec theta, arma::vec apars, double linpred1, double linpred2, double dpnorm1, double dpnorm2, arma::uword mj, arma::uword pj, arma::uword nparj){
    //  SJ: Compute d2log(Pr) / dtheta dpar^T
    arma::cube d2logPidtdu = zeros<cube>(pj, nparj, 1);
    for(arma::uword uu = 0; uu < pj; uu++){
        //  SJ: a-par
        for(arma::uword temppar = 0; temppar < pj; temppar++){
            d2logPidtdu(uu, temppar, 0) = (linpred1 * dpnorm1 - linpred2 * dpnorm2 + pow(dpnorm1 - dpnorm2,2.0)) * apars(uu) * theta(temppar);
        }
        
        //  SJ: b-par
        if(ynj == 1){
            d2logPidtdu(uu, pj, 0) = -(linpred2 * dpnorm2 + pow(dpnorm2,2.0)) * apars(uu); 
        }
        else if(ynj == mj){
            d2logPidtdu(uu, ynj - 2 + pj, 0) = (linpred1 * dpnorm1 - pow(dpnorm1,2.0)) * apars(uu); 
        }
        else{
            d2logPidtdu(uu, ynj - 2 + pj, 0) = (linpred1 * dpnorm1 - pow(dpnorm1,2.0) + dpnorm1 * dpnorm2) * apars(uu); 
            d2logPidtdu(uu, ynj - 1 + pj, 0) = (-linpred2 * dpnorm2 - pow(dpnorm2,2.0) + dpnorm1 * dpnorm2) * apars(uu); 
        }

    }
    
    return(d2logPidtdu);
    
}

// [[Rcpp::export]]
double d3gjd3t_GRM_probit(arma::uword ynj, arma::vec apars, double linpred1, double linpred2, double dpnorm1, double dpnorm2, arma::uword mj, arma::uword pj){
    //  SJ: Compute the constant part in d2logPi / d2theta.
    //      In order to get d2logPi / d2theta, we multiply the output by corresponding entries from apar.
    double d3logPid3t = 0.0;
    if(ynj == 1){
        d3logPid3t = dpnorm2 - pow(linpred2,2.0) * dpnorm2 - 3.0 * linpred2 * pow(dpnorm2,2.0) - 2.0 * pow(dpnorm2,3.0) ;
    }
    else if(ynj == mj){
        d3logPid3t = -1.0*dpnorm1 + pow(linpred1,2.0) * dpnorm1 - 3.0 * linpred1 * pow(dpnorm1,2.0) + 2.0 * pow(dpnorm1,3.0) ;
    }
    else {
        d3logPid3t = -(dpnorm1 - dpnorm2) - pow(linpred2,2.0) * dpnorm2 + pow(linpred1,2.0) * dpnorm1 - 3.0 * (dpnorm1 - dpnorm2) * ( linpred1 * dpnorm1 - linpred2 * dpnorm2 ) + 2.0*pow(dpnorm1 - dpnorm2,3); 
    }
    return(d3logPid3t);
}

// [[Rcpp::export]]
double d4gjd4t_GRM_probit(arma::uword ynj, arma::vec apars, double linpred1, double linpred2, double dpnorm1, double dpnorm2, arma::uword mj, arma::uword pj){
    //  SJ: Compute the constant part in d2logPi / d2theta.
    //      In order to get d2logPi / d2theta, we multiply the output by corresponding entries from apar.
    double d4logPid4t = 0.0;
    if(ynj == 1){
        d4logPid4t = 3.0 * linpred2 * dpnorm2 - pow(linpred2,3)*dpnorm2 - 6.0 * pow(dpnorm2,4) - 3.0 * pow(linpred2 * dpnorm2,2) + 4.0*( dpnorm2 - pow(linpred2,2)*dpnorm2 )*dpnorm2  - 12.0 * linpred2 * pow(dpnorm2,3 ) ;
    }
    else if(ynj == mj){
        d4logPid4t = -3.0 * linpred1 * dpnorm1 + pow(linpred1,3)*dpnorm1 - 6.0 * pow(dpnorm1,4) - 3.0 * pow(linpred1 * dpnorm1,2) - 4.0*( -1.0 * dpnorm1 + pow(linpred1,2)*dpnorm1 )*dpnorm1 + 12.0*( linpred1 * dpnorm1 )*pow(dpnorm1,2 ) ;
    }
    else {
        // SJ: This formula is a general expression, even valid for 1 and mj.
        d4logPid4t = 3.0 * linpred2 * dpnorm2 - 3.0 * linpred1 * dpnorm1 - 3.0 * pow( linpred1 * dpnorm1 - linpred2 * dpnorm2, 2.0 ) - 4.0*pow(linpred1 * dpnorm1, 2.0) - 4.0*pow(linpred2 * dpnorm2, 2.0) - pow(linpred2, 3.0) * dpnorm2 + pow(linpred1, 3.0) * dpnorm1 - 6.0 * pow(dpnorm1 - dpnorm2 , 4.0) + 4.0 * pow(dpnorm2, 2.0) + 4.0 * pow(dpnorm1, 2.0) + 12.0*( linpred1 * dpnorm1 - linpred2 * dpnorm2 )*pow(dpnorm1 - dpnorm2, 2.0) + 4.0 * pow(linpred1,2) * dpnorm1 * dpnorm2 + 4.0 * pow(linpred2,2) * dpnorm1 * dpnorm2 - 8.0 * dpnorm1 * dpnorm2;
    }
    return(d4logPid4t);
}

// [[Rcpp::export]]
double d5gjd5t_GRM_probit(arma::uword ynj, arma::vec apars, double linpred1, double linpred2, double dpnorm1, double dpnorm2, arma::uword mj, arma::uword pj){
    //  SJ: Compute the constant part in d2logPi / d2theta.
    //      In order to get d2logPi / d2theta, we multiply the output by corresponding entries from apar.
    double d5logPid5t = 0.0;
    if(ynj == 1){
        d5logPid5t = -3.0 * dpnorm2 + 6.0 * pow(linpred2, 2.0) * dpnorm2 + 25.0 * linpred2 * pow(dpnorm2, 2.0) + 20.0 * pow(dpnorm2, 3.0) - pow(linpred2, 4.0) * dpnorm2 - 15.0 * pow(linpred2, 3.0) * pow(dpnorm2, 2.0) - 50.0 * pow(linpred2, 2.0) * pow(dpnorm2, 3.0) - 60.0 * linpred2 * pow(dpnorm2, 4.0) - 24.0 * pow(dpnorm2, 5.0) ;
        }
    else if(ynj == mj){
        d5logPid5t = 3.0 * dpnorm1 - 6.0 * pow(linpred1, 2.0) * dpnorm1 + 25.0 * linpred1 * pow(dpnorm1, 2.0) - 20.0 * pow(dpnorm1, 3.0) + pow(linpred1, 4.0) * dpnorm1 - 15.0 * pow(linpred1, 3.0) * pow(dpnorm1, 2.0) + 50.0 * pow(linpred1, 2.0) * pow(dpnorm1, 3.0) - 60.0 * linpred1 * pow(dpnorm1, 4.0) + 24.0 * pow(dpnorm1, 5.0) ;
    }
    else {
        // SJ: This formula is a general expression, even valid for 1 and mj.
        d5logPid5t = 3.0 * (dpnorm1 - dpnorm2) - 6.0 * pow(linpred1,2) * dpnorm1 + 6.0 * pow(linpred2,2) * dpnorm2 + pow(linpred1,4) * dpnorm1 - pow(linpred2,4) * dpnorm2 - 5.0 * ( 3.0 * linpred2 - pow(linpred2, 3.0) ) * (dpnorm1 - dpnorm2) * dpnorm2 + 5.0 * ( 3.0 * linpred1 - pow(linpred1, 3.0) ) * (dpnorm1 - dpnorm2) * dpnorm1 + 10.0 * (linpred1 * dpnorm1 - linpred2 * dpnorm2) * ( dpnorm1 - pow(linpred1, 2.0) * dpnorm1 - dpnorm2 + pow(linpred2, 2.0) * dpnorm2) + 30.0 * pow(linpred1 * dpnorm1 - linpred2 * dpnorm2, 2.0) * (dpnorm1 - dpnorm2) - 20.0 * pow(dpnorm1 - dpnorm2, 3.0) - 20.0 * pow(linpred2, 2.0) * dpnorm2 * pow(dpnorm1 - dpnorm2, 2.0) + 20.0 * pow(linpred1, 2.0) * dpnorm1 * pow(dpnorm1 - dpnorm2, 2.0) - 60.0 * pow(dpnorm1 - dpnorm2, 3.0) * (linpred1 * dpnorm1 - linpred2 * dpnorm2) + 24.0 * pow(dpnorm1 - dpnorm2, 5.0) ;
    }
    
    return(d5logPid5t);
}

// [[Rcpp::export]]
Rcpp::List item_GRM( arma::uword ynj, arma::vec thetajj, arma::vec aparsjj, arma::vec bparsjj, std::string link, arma::uword mj_jj, arma::uword pj, arma::uword nparj ){
    // SJ: Using Rcpp::List can be slow. I am using mat instead, but the indices are not straightforward. 
    //We only need two cumulative probabilities for the GRM and, thus, only the derivatives of these two.
    //mj_jj is the number of response categories.
    //We need the entries that correspond to ynj and ynj + 1.
	double linpred1 = 0.0;
	double linpred2 = 0.0;
	double exp1 = 0.0;
	double exp2 = 0.0;
	double pdfnorm1 = 0.0;
	double pdfnorm2 = 0.0;
	double cdfnorm1 = 0.0;
	double cdfnorm2 = 0.0;
	double sumat;
    if(link == "logit"){
        
        arma::vec Pi = zeros<vec>(2);
        arma::mat dPidt = zeros<mat>(pj, 2);
        arma::cube d2Pidt2 = zeros<cube>(pj, pj, 2);
        arma::mat dPidu = zeros<mat>(nparj, 2);
        arma::cube d2Pidtdu = zeros<cube>(pj, nparj, 2);
        sumat = sum(aparsjj % thetajj);
        arma::mat dlogPrdt = zeros<mat>(pj,1) ;
        arma::mat d2logPrdt2 = zeros<mat>(pj,pj) ;
        arma::cube d3logPrdt3 = zeros<cube>(pj,pj,pj) ;
        
        if(ynj == 1){
            linpred2 = sumat + bparsjj(ynj);
            exp2 = exp(-linpred2);
            Pi(0) = 1.0;
            Pi(1) = 1.0 / (1.0 + exp2);
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 1) = pow(Pi(1), 2) * exp2 * thetajj(temppar);
            //b-par
            dPidu(pj, 1) = pow(Pi(1), 2) * exp2;
        }
        else if(ynj == mj_jj){
            linpred1 = sumat + bparsjj(ynj - 1);
            exp1 = exp(-linpred1);
            Pi(0) = 1.0 / (1.0 + exp1);
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 0) = pow(Pi(0), 2) * exp1 * thetajj(temppar);
            //b-par
            dPidu(ynj - 2 + pj, 0) = pow(Pi(0), 2) * exp1;
        }
        else{
            linpred1 = sumat + bparsjj(ynj - 1);
            linpred2 = sumat + bparsjj(ynj);
            exp1 = exp(-linpred1);
            exp2 = exp(-linpred2);
            Pi(0) = 1.0 / (1.0 + exp1);
            Pi(1) = 1.0 / (1.0 + exp2);
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++){
                dPidu(temppar, 0) = pow(Pi(0), 2) * exp1 * thetajj(temppar);
                dPidu(temppar, 1) = pow(Pi(1), 2) * exp2 * thetajj(temppar);
            }
            //b-par
            dPidu(ynj - 2 + pj, 0) = pow(Pi(0), 2) * exp1;
            dPidu(ynj - 1 + pj, 1) = pow(Pi(1), 2) * exp2;
        }
        
        // SJ: Gradient of log(Prob) w.r.t. theta
        for(arma::uword tt = 0; tt < 2; tt++){
            for(arma::uword uu = 0; uu < pj; uu++){
                dPidt(uu, tt) = aparsjj(uu) * Pi(tt) * (1.0 - Pi(tt));
            }
            for(arma::uword uu = 0; uu < pj; uu++){
                for(arma::uword temppar = 0; temppar < nparj; temppar++){
                    d2Pidtdu(uu, temppar, tt) = aparsjj(uu) * (dPidu(temppar, tt) * (1.0 - Pi(tt)) - Pi(tt) * dPidu(temppar, tt));
                }
                d2Pidtdu(uu, uu, tt) += Pi(tt) * (1.0 - Pi(tt));
                for(arma::uword vv = 0; vv < pj; vv++){
                    d2Pidt2(uu, vv, tt) = aparsjj(uu) * dPidt(vv, tt) * (1.0 - 2.0 * Pi(tt));
                }
            }
        }
        
        double Pr = Pi(0) - Pi(1);
        double logPr = log(Pr);
        double d1hconst = 1.0 - Pi(0) - Pi(1);
        double d2hconst0 = Pi(0) * (1.0 - Pi(0));
        double d2hconst1 = Pi(1) * (1.0 - Pi(1));
        double d3hconst0 = Pi(0) * (1.0 - Pi(0)) * (1.0 - 2.0*Pi(0));
        double d3hconst1 = Pi(1) * (1.0 - Pi(1)) * (1.0 - 2.0*Pi(1));
        for(arma::uword uu = 0; uu < pj; uu++){
            dlogPrdt(uu, 0) = aparsjj(uu) * d1hconst;
            for(arma::uword vv = 0; vv < pj; vv++){
                d2logPrdt2(uu, vv) = -1.0 * aparsjj(uu) * aparsjj(vv) * (d2hconst0 + d2hconst1);
                for(arma::uword ww = 0; ww < pj; ww++){
                    d3logPrdt3(uu, vv, ww) = -1.0 * aparsjj(uu) * aparsjj(vv) * aparsjj(ww) * (d3hconst0 + d3hconst1);
                }
            }
        }
        
        //  Gradient Gradient of log(Prob) w.r.t. item parameters
        arma::mat dlogPrdu = zeros<mat>(nparj,1) ;
        arma::mat d2logPrdtdu = zeros<mat>(pj, nparj);
        arma::cube d3logPrdt2du = zeros<cube>(pj,pj,nparj) ;
        //a-par
        for(arma::uword temppar = 0; temppar < pj; temppar++) dlogPrdu(temppar, 0) = d1hconst * thetajj(temppar);
        //  Cross derivative of log(Prob) w.r.t. latent variable and apar
        for(arma::uword uu = 0; uu < pj; uu++){
            for(arma::uword temppar = 0; temppar < pj; temppar++){
                d2logPrdtdu(uu, temppar) = -1.0 * aparsjj(uu) * thetajj(temppar) * (d2hconst0 + d2hconst1);
                if( uu == temppar ) d2logPrdtdu(uu, temppar) += d1hconst;
            }
            
            for(arma::uword vv = 0; vv < pj; vv++){
                for(arma::uword temppar = 0; temppar < pj; temppar++){
                    d3logPrdt2du(uu, vv, temppar) = -1.0 * aparsjj(uu) * aparsjj(vv) * thetajj(temppar) * (d3hconst0 + d3hconst1);
                    if( uu == temppar ) d3logPrdt2du(uu, vv, temppar) += -1.0 * aparsjj(vv) * (d2hconst0 + d2hconst1);
                    if( vv == temppar ) d3logPrdt2du(uu, vv, temppar) += -1.0 * aparsjj(uu) * (d2hconst0 + d2hconst1);
                }
            }
        }

        //  bpar
        if(ynj == 1){
            dlogPrdu(pj, 0) = d1hconst; //-1.0 * Pi(1);
            for(arma::uword uu = 0; uu < pj; uu++){
                d2logPrdtdu(uu, pj) = -1.0 * aparsjj(uu) * d2hconst1; // SJ: I think we should have -1 here.
                for(arma::uword vv = 0; vv < pj; vv++){
                    d3logPrdt2du(uu, vv, pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst1;
                }
            }
        }
        else if(ynj == mj_jj){
            dlogPrdu(ynj - 2 + pj, 0) = d1hconst; //1.0 - Pi(0);
            for(arma::uword uu = 0; uu < pj; uu++){
                d2logPrdtdu(uu, ynj - 2 + pj) = -1.0 * aparsjj(uu) * d2hconst0; //  SJ: I think we should have -1 here.
                for(arma::uword vv = 0; vv < pj; vv++){
                    d3logPrdt2du(uu, vv, ynj - 2 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst0;
                }
            }
        }
        else{
            dlogPrdu(ynj - 2 + pj, 0) = d2hconst0 / (Pi(0) - Pi(1));  
            dlogPrdu(ynj - 1 + pj, 0) = -1.0 * d2hconst1 / (Pi(0) - Pi(1));
            for(arma::uword uu = 0; uu < pj; uu++){
                d2logPrdtdu(uu, ynj - 2 + pj) = -1.0 * aparsjj(uu) * d2hconst0; // SJ: I think we should have -1 here.
                d2logPrdtdu(uu, ynj - 1 + pj) = -1.0 * aparsjj(uu) * d2hconst1; // SJ: I think we should have -1 here.
                for(arma::uword vv = 0; vv < pj; vv++){
                    d3logPrdt2du(uu, vv, ynj - 2 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst0;
                    d3logPrdt2du(uu, vv, ynj - 1 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst1;
                }
            }
        }
        
        return Rcpp::List::create(Rcpp::Named("Pi") = Pi,
                                  Rcpp::Named("dPidt") = dPidt,
                                  Rcpp::Named("d2Pidt2") = d2Pidt2,
                                  Rcpp::Named("dPidu") = dPidu,
                                  Rcpp::Named("d2Pidtdu") = d2Pidtdu,
                                  Rcpp::Named("sumat") = sumat,
                                  Rcpp::Named("linpred1") = linpred1,
                                  Rcpp::Named("linpred2") = linpred2,
                                  Rcpp::Named("exp1") = exp1,
                                  Rcpp::Named("exp2") = exp2,
                                  Rcpp::Named("logPr") = logPr,
                                  Rcpp::Named("dlogPrdt") = dlogPrdt,
                                  Rcpp::Named("d2logPrdt2") = d2logPrdt2,
                                  Rcpp::Named("d3logPrdt3") = d3logPrdt3,
                                  Rcpp::Named("dlogPrdu") = dlogPrdu,
                                  Rcpp::Named("d2logPrdtdu") = d2logPrdtdu,
                                  Rcpp::Named("d3logPrdt2du") = d3logPrdt2du);
        
    } 
    else {
        //  SJ: For the probit link, we need both dnorm() and pnorm().
        arma::vec Pi = zeros<vec>(2);
        arma::mat dPidt = zeros<mat>(pj, 2);
        arma::cube d2Pidt2 = zeros<cube>(pj, pj, 2);
        arma::mat dPidu = zeros<mat>(nparj, 2);
        arma::cube d2Pidtdu = zeros<cube>(pj, nparj, 2);
		sumat = sum(aparsjj % thetajj);
        
        
        if(ynj == 1){
            linpred2 = sumat + bparsjj(ynj);
            pdfnorm2 = R::dnorm( -linpred2, 0.0, 1.0, 0L ); 
            cdfnorm2 = R::pnorm( -linpred2, 0.0, 1.0, 0L, 0L ); //  lower=false.
            Pi(0) = 1.0;
            Pi(1) = cdfnorm2;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 1) = pdfnorm2 * thetajj(temppar);
            //b-par
            dPidu(pj, 1) = pdfnorm2;
        }
        else if(ynj == mj_jj){
            linpred1 = sumat + bparsjj(ynj - 1);
            pdfnorm1 = R::dnorm( -linpred1, 0.0, 1.0, 0L ); 
            cdfnorm1 = R::pnorm( -linpred1, 0.0, 1.0, 0L, 0L );
            Pi(0) = cdfnorm1;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 0) = pdfnorm1 * thetajj(temppar);
            //b-par
            dPidu(ynj - 2 + pj, 0) = pdfnorm1;
        }
        else{
            linpred1 = sumat + bparsjj(ynj - 1);
            linpred2 = sumat + bparsjj(ynj);
            pdfnorm1 = R::dnorm( -linpred1, 0.0, 1.0, 0L ); 
            cdfnorm1 = R::pnorm( -linpred1, 0.0, 1.0, 0L, 0L );
            pdfnorm2 = R::dnorm( -linpred2, 0.0, 1.0, 0L ); 
            cdfnorm2 = R::pnorm( -linpred2, 0.0, 1.0, 0L, 0L );
            Pi(0) = cdfnorm1;
            Pi(1) = cdfnorm2;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++){
                dPidu(temppar, 0) = pdfnorm1 * thetajj(temppar);
                dPidu(temppar, 1) = pdfnorm2 * thetajj(temppar);
            }
            //b-par
            dPidu(ynj - 2 + pj, 0) = pdfnorm1;
            dPidu(ynj - 1 + pj, 1) = pdfnorm2;
        }
        
        //  SJ: Compute dPi / dtheta and d2Pi / d2theta
        for(arma::uword tt = 0; tt < 2; tt++){
            for(arma::uword uu = 0; uu < pj; uu++){
                dPidt(uu, tt) = (pdfnorm1 - pdfnorm2) / (Pi(0) - Pi(1)) * aparsjj(uu) ;
                for(arma::uword temppar = 0; temppar < nparj; temppar++){
                    d2Pidtdu(uu, temppar, tt) = aparsjj(uu) * (dPidu(temppar, tt) * (1.0 - Pi(tt)) - Pi(tt) * dPidu(temppar, tt));
                }
                d2Pidtdu(uu, uu, tt) += Pi(tt) * (1.0 - Pi(tt));
                for(arma::uword vv = 0; vv < pj; vv++){
                    d2Pidt2(uu, vv, tt) = aparsjj(uu) * dPidt(vv, tt) * (1.0 - 2.0 * Pi(tt));
                }
            }
        }
        
        return Rcpp::List::create(Rcpp::Named("Pi") = Pi,
                                  Rcpp::Named("dPidt") = dPidt,
                                  Rcpp::Named("d2Pidt2") = d2Pidt2,
                                  Rcpp::Named("dPidu") = dPidu,
                                  Rcpp::Named("d2Pidtdu") = d2Pidtdu,
                                  Rcpp::Named("sumat") = sumat,
                                  Rcpp::Named("linpred1") = linpred1,
                                  Rcpp::Named("linpred2") = linpred2,
                                  Rcpp::Named("exp1") = exp1,
                                  Rcpp::Named("exp2") = exp2);
        
    }
    
    
    
}

// [[Rcpp::export]]
arma::mat item_GRM_quad( arma::uword ynj, arma::vec thetajj, arma::vec aparsjj, arma::vec bparsjj, std::string link, arma::uword mj_jj, arma::uword pj, arma::uword nparj ){
    // SJ: Using Rcpp::List can be slow. I am using mat instead, but the indices are not straightforward. 
    //We only need two cumulative probabilities for the GRM and, thus, only the derivatives of these two.
    //mj_jj is the number of response categories.
    //We need the entries that correspond to ynj and ynj + 1.
	double linpred1;
	double linpred2;
	//double exp1;
	//double exp2;
	double pdfnorm1 = 0.0;
	double pdfnorm2 = 0.0;
	double cdfnorm1 = 0.0;
	double cdfnorm2 = 0.0;
	double sumat;
    if(link == "logit"){
        
        // SJ: In order to avoid using lists, we will return a matrix. The matrix will be organized as
        //     | Pi(0)  dlog(Pr) / dtheta  d2log(Pr) / dtheta^2
        //     | Pi(1)         ...                ...
        //     |  0     dlog(Pr) / dtheta  d2log(Pr) / dtheta^2
        //     |  0     dlog(Pr) / du      d2log(Pr) / dudtheta
        //     |  ...         ...                  ...
        //     |  0     dlog(Pr) / du      d2log(Pr) / dudtheta
        //     such that Pr = Pi(0) - Pi(1)
        arma::mat ret;
        if( (pj + nparj) > 2 ){
            ret = zeros<mat>(pj + nparj, 1 + 1 + pj) ;
        } else {
            ret = zeros<mat>(2, 1 + 1 + pj) ;
        }
        arma::vec Pi = zeros<vec>(2);
        arma::mat dPidt = zeros<mat>(pj, 2);
        arma::cube d2Pidt2 = zeros<cube>(pj, pj, 2);
        arma::mat dPidu = zeros<mat>(nparj, 2);
        arma::cube d2Pidtdu = zeros<cube>(pj, nparj, 2);
        double sumat = sum(aparsjj % thetajj);
        double linpred1;
        double linpred2;
        double exp1;
        double exp2;
        arma::mat dlogPrdt = zeros<mat>(pj,1) ;
        arma::mat d2logPrdt2 = zeros<mat>(pj,pj) ;
        arma::cube d3logPrdt3 = zeros<cube>(pj,pj,pj) ;
        
        if(ynj == 1){
            linpred2 = sumat + bparsjj(ynj);
            exp2 = exp(-linpred2);
            Pi(0) = 1.0;
            Pi(1) = 1.0 / (1.0 + exp2);
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 1) = pow(Pi(1), 2) * exp2 * thetajj(temppar);
            //b-par
            dPidu(pj, 1) = pow(Pi(1), 2) * exp2;
        }
        else if(ynj == mj_jj){
            linpred1 = sumat + bparsjj(ynj - 1);
            exp1 = exp(-linpred1);
            Pi(0) = 1.0 / (1.0 + exp1);
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 0) = pow(Pi(0), 2) * exp1 * thetajj(temppar);
            //b-par
            dPidu(ynj - 2 + pj, 0) = pow(Pi(0), 2) * exp1;
        }
        else{
            linpred1 = sumat + bparsjj(ynj - 1);
            linpred2 = sumat + bparsjj(ynj);
            exp1 = exp(-linpred1);
            exp2 = exp(-linpred2);
            Pi(0) = 1.0 / (1.0 + exp1);
            Pi(1) = 1.0 / (1.0 + exp2);
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++){
                dPidu(temppar, 0) = pow(Pi(0), 2) * exp1 * thetajj(temppar);
                dPidu(temppar, 1) = pow(Pi(1), 2) * exp2 * thetajj(temppar);
            }
            //b-par
            dPidu(ynj - 2 + pj, 0) = pow(Pi(0), 2) * exp1;
            dPidu(ynj - 1 + pj, 1) = pow(Pi(1), 2) * exp2;
        }
        
        // SJ: Gradient of log(Prob) w.r.t. theta
        for(arma::uword tt = 0; tt < 2; tt++){
            for(arma::uword uu = 0; uu < pj; uu++){
                dPidt(uu, tt) = aparsjj(uu) * Pi(tt) * (1.0 - Pi(tt));
            }
            for(arma::uword uu = 0; uu < pj; uu++){
                for(arma::uword temppar = 0; temppar < nparj; temppar++){
                    d2Pidtdu(uu, temppar, tt) = aparsjj(uu) * (dPidu(temppar, tt) * (1.0 - Pi(tt)) - Pi(tt) * dPidu(temppar, tt));
                }
                d2Pidtdu(uu, uu, tt) += Pi(tt) * (1.0 - Pi(tt));
                for(arma::uword vv = 0; vv < pj; vv++){
                    d2Pidt2(uu, vv, tt) = aparsjj(uu) * dPidt(vv, tt) * (1.0 - 2.0 * Pi(tt));
                }
            }
        }
        
        //double Pr = Pi(0) - Pi(1);
        //double logPr = log(Pr);
        double d1hconst = 1.0 - Pi(0) - Pi(1);
        double d2hconst0 = Pi(0) * (1.0 - Pi(0));
        double d2hconst1 = Pi(1) * (1.0 - Pi(1));
        double d3hconst0 = Pi(0) * (1.0 - Pi(0)) * (1.0 - 2.0*Pi(0));
        double d3hconst1 = Pi(1) * (1.0 - Pi(1)) * (1.0 - 2.0*Pi(1));
        for(arma::uword uu = 0; uu < pj; uu++){
            dlogPrdt(uu, 0) = aparsjj(uu) * d1hconst;
            for(arma::uword vv = 0; vv < pj; vv++){
                d2logPrdt2(uu, vv) = -1.0 * aparsjj(uu) * aparsjj(vv) * (d2hconst0 + d2hconst1);
                for(arma::uword ww = 0; ww < pj; ww++){
                    d3logPrdt3(uu, vv, ww) = -1.0 * aparsjj(uu) * aparsjj(vv) * aparsjj(ww) * (d3hconst0 + d3hconst1);
                }
            }
        }
        
        //  Gradient Gradient of log(Prob) w.r.t. item parameters
        arma::mat dlogPrdu = zeros<mat>(nparj,1) ;
        arma::mat d2logPrdtdu = zeros<mat>(pj, nparj);
        arma::cube d3logPrdt2du = zeros<cube>(pj,pj,nparj) ;
        //a-par
        for(arma::uword temppar = 0; temppar < pj; temppar++) dlogPrdu(temppar, 0) = d1hconst * thetajj(temppar);
        //  Cross derivative of log(Prob) w.r.t. latent variable and apar
        for(arma::uword uu = 0; uu < pj; uu++){
            for(arma::uword temppar = 0; temppar < pj; temppar++){
                d2logPrdtdu(uu, temppar) = -1.0 * aparsjj(uu) * thetajj(temppar) * (d2hconst0 + d2hconst1);
                if( uu == temppar ) d2logPrdtdu(uu, temppar) += d1hconst;
            }
            
            for(arma::uword vv = 0; vv < pj; vv++){
                for(arma::uword temppar = 0; temppar < pj; temppar++){
                    d3logPrdt2du(uu, vv, temppar) = -1.0 * aparsjj(uu) * aparsjj(vv) * thetajj(temppar) * (d3hconst0 + d3hconst1);
                    if( uu == temppar ) d3logPrdt2du(uu, vv, temppar) += -1.0 * aparsjj(vv) * (d2hconst0 + d2hconst1);
                    if( vv == temppar ) d3logPrdt2du(uu, vv, temppar) += -1.0 * aparsjj(uu) * (d2hconst0 + d2hconst1);
                }
            }
        }
        
        //  bpar
        if(ynj == 1){
            dlogPrdu(pj, 0) = d1hconst; //-1.0 * Pi(1);
            for(arma::uword uu = 0; uu < pj; uu++){
                d2logPrdtdu(uu, pj) = -1.0 * aparsjj(uu) * d2hconst1; // SJ: I think we should have -1 here.
                for(arma::uword vv = 0; vv < pj; vv++){
                    d3logPrdt2du(uu, vv, pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst1;
                }
            }
        }
        else if(ynj == mj_jj){
            dlogPrdu(ynj - 2 + pj, 0) = d1hconst; //1.0 - Pi(0);
            for(arma::uword uu = 0; uu < pj; uu++){
                d2logPrdtdu(uu, ynj - 2 + pj) = -1.0 * aparsjj(uu) * d2hconst0; //  SJ: I think we should have -1 here.
                for(arma::uword vv = 0; vv < pj; vv++){
                    d3logPrdt2du(uu, vv, ynj - 2 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst0;
                }
            }
        }
        else{
            dlogPrdu(ynj - 2 + pj, 0) = d2hconst0 / (Pi(0) - Pi(1));  
            dlogPrdu(ynj - 1 + pj, 0) = -1.0 * d2hconst1 / (Pi(0) - Pi(1));
            for(arma::uword uu = 0; uu < pj; uu++){
                d2logPrdtdu(uu, ynj - 2 + pj) = -1.0 * aparsjj(uu) * d2hconst0; // SJ: I think we should have -1 here.
                d2logPrdtdu(uu, ynj - 1 + pj) = -1.0 * aparsjj(uu) * d2hconst1; // SJ: I think we should have -1 here.
                for(arma::uword vv = 0; vv < pj; vv++){
                    d3logPrdt2du(uu, vv, ynj - 2 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst0;
                    d3logPrdt2du(uu, vv, ynj - 1 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst1;
                }
            }
        }
        
        ret(0, 0) = Pi(0);
        ret(1, 0) = Pi(1);
        ret.submat(0, 1, pj + nparj -1, 1) = join_cols(dlogPrdt, dlogPrdu);
        ret.submat(0, 2, pj + nparj -1, 1 + pj) = join_cols(d2logPrdt2, d2logPrdtdu.t());
        
        return(ret);
        
    } 
    else {
        //  SJ: For the probit link, we need both dnorm() and pnorm().
        arma::vec Pi = zeros<vec>(2);
        arma::mat dPidt = zeros<mat>(pj, 2);
        arma::cube d2Pidt2 = zeros<cube>(pj, pj, 2);
        arma::mat dPidu = zeros<mat>(nparj, 2);
        arma::cube d2Pidtdu = zeros<cube>(pj, nparj, 2);
        sumat = sum(aparsjj % thetajj);

        
        if(ynj == 1){
            linpred2 = sumat + bparsjj(ynj);
            pdfnorm2 = R::dnorm( -linpred2, 0.0, 1.0, 0L ); 
            cdfnorm2 = R::pnorm( -linpred2, 0.0, 1.0, 0L, 0L ); //  lower=false.
            Pi(0) = 1.0;
            Pi(1) = cdfnorm2;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 1) = pdfnorm2 * thetajj(temppar);
            //b-par
            dPidu(pj, 1) = pdfnorm2;
        }
        else if(ynj == mj_jj){
            linpred1 = sumat + bparsjj(ynj - 1);
            pdfnorm1 = R::dnorm( -linpred1, 0.0, 1.0, 0L ); 
            cdfnorm1 = R::pnorm( -linpred1, 0.0, 1.0, 0L, 0L );
            Pi(0) = cdfnorm1;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 0) = pdfnorm1 * thetajj(temppar);
            //b-par
            dPidu(ynj - 2 + pj, 0) = pdfnorm1;
        }
        else{
            linpred1 = sumat + bparsjj(ynj - 1);
            linpred2 = sumat + bparsjj(ynj);
            pdfnorm1 = R::dnorm( -linpred1, 0.0, 1.0, 0L ); 
            cdfnorm1 = R::pnorm( -linpred1, 0.0, 1.0, 0L, 0L );
            pdfnorm2 = R::dnorm( -linpred2, 0.0, 1.0, 0L ); 
            cdfnorm2 = R::pnorm( -linpred2, 0.0, 1.0, 0L, 0L );
            Pi(0) = cdfnorm1;
            Pi(1) = cdfnorm2;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++){
                dPidu(temppar, 0) = pdfnorm1 * thetajj(temppar);
                dPidu(temppar, 1) = pdfnorm2 * thetajj(temppar);
            }
            //b-par
            dPidu(ynj - 2 + pj, 0) = pdfnorm1;
            dPidu(ynj - 1 + pj, 1) = pdfnorm2;
        }
        
        //  SJ: Compute dPi / dtheta and d2Pi / d2theta
        for(arma::uword tt = 0; tt < 2; tt++){
            for(arma::uword uu = 0; uu < pj; uu++){
                dPidt(uu, tt) = (pdfnorm1 - pdfnorm2) / (Pi(0) - Pi(1)) * aparsjj(uu) ;
                for(arma::uword temppar = 0; temppar < nparj; temppar++){
                    d2Pidtdu(uu, temppar, tt) = aparsjj(uu) * (dPidu(temppar, tt) * (1.0 - Pi(tt)) - Pi(tt) * dPidu(temppar, tt));
                }
                d2Pidtdu(uu, uu, tt) += Pi(tt) * (1.0 - Pi(tt));
                for(arma::uword vv = 0; vv < pj; vv++){
                    d2Pidt2(uu, vv, tt) = aparsjj(uu) * dPidt(vv, tt) * (1.0 - 2.0 * Pi(tt));
                }
            }
        }
        
        return(dPidu);
        
    }
    
}

// [[Rcpp::export]]
arma::cube item_GRM_3rd( arma::uword ynj, arma::vec thetajj, arma::vec aparsjj, arma::vec bparsjj, std::string link, arma::uword mj_jj, arma::uword pj, arma::uword nparj, arma::vec Pi, bool dt2du ){
    // SJ: This function returns d3logPrdt3, possibly also d3logPrdt2du, which is controlled by dt2du.
    //     If dt2du == true, then we also return d3logPrdt2du.
    //We only need two cumulative probabilities for the GRM and, thus, only the derivatives of these two.
    //mj_jj is the number of response categories.
    //We need the entries that correspond to ynj and ynj + 1.
    arma::cube ret;
	double linpred1;
	double linpred2;
	//double exp1;
	//double exp2;
	double pdfnorm1 = 0.0;
	double pdfnorm2 = 0.0;
	double cdfnorm1 = 0.0;
	double cdfnorm2 = 0.0;
	double sumat;
    if( dt2du == true ){
        ret = zeros<cube>(pj, pj, pj + nparj) ;
    } else {
        ret = zeros<cube>(pj, pj, pj) ;
    }
    if(link == "logit"){

        arma::cube d3logPrdt3 = zeros<cube>(pj,pj,pj) ;
        double d2hconst0 = Pi(0) * (1.0 - Pi(0));
        double d2hconst1 = Pi(1) * (1.0 - Pi(1));
        double d3hconst0 = Pi(0) * (1.0 - Pi(0)) * (1.0 - 2.0*Pi(0));
        double d3hconst1 = Pi(1) * (1.0 - Pi(1)) * (1.0 - 2.0*Pi(1));
        for(arma::uword uu = 0; uu < pj; uu++){
            for(arma::uword vv = 0; vv < pj; vv++){
                for(arma::uword ww = 0; ww < pj; ww++){
                    d3logPrdt3(uu, vv, ww) = -1.0 * aparsjj(uu) * aparsjj(vv) * aparsjj(ww) * (d3hconst0 + d3hconst1);
                }
            }
        }

        arma::cube d3logPrdt2du = zeros<cube>(pj,pj,nparj) ;
        for(arma::uword uu = 0; uu < pj; uu++){
            for(arma::uword vv = 0; vv < pj; vv++){
                // SJ: apar
                for(arma::uword temppar = 0; temppar < pj; temppar++){
                    d3logPrdt2du(uu, vv, temppar) = -1.0 * aparsjj(uu) * aparsjj(vv) * thetajj(temppar) * (d3hconst0 + d3hconst1);
                    if( uu == temppar ) d3logPrdt2du(uu, vv, temppar) += -1.0 * aparsjj(vv) * (d2hconst0 + d2hconst1);
                    if( vv == temppar ) d3logPrdt2du(uu, vv, temppar) += -1.0 * aparsjj(uu) * (d2hconst0 + d2hconst1);
                }
                // SJ: bpar
                if(ynj == 1){
                    d3logPrdt2du(uu, vv, pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst1;
                }
                else if(ynj == mj_jj){
                    d3logPrdt2du(uu, vv, ynj - 2 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst0;
                }
                else{
                    d3logPrdt2du(uu, vv, ynj - 2 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst0;
                    d3logPrdt2du(uu, vv, ynj - 1 + pj) = -1.0 * aparsjj(uu) * aparsjj(vv) * d3hconst1;
                }
            }
        }
        
        ret( span(0, pj - 1), span(0, pj - 1), span(0, pj - 1) ) = d3logPrdt3;
        if( dt2du == true ){
            ret( span(0, pj - 1), span(0, pj - 1), span(pj, pj + nparj - 1) ) = d3logPrdt2du;
        } 
        
    } 
    else {
        //  SJ: For the probit link, we need both dnorm() and pnorm().
        arma::vec Pi = zeros<vec>(2);
        arma::mat dPidt = zeros<mat>(pj, 2);
        arma::cube d2Pidt2 = zeros<cube>(pj, pj, 2);
        arma::mat dPidu = zeros<mat>(nparj, 2);
        arma::cube d2Pidtdu = zeros<cube>(pj, nparj, 2);
        sumat = sum(aparsjj % thetajj);
       
        
        if(ynj == 1){
            linpred2 = sumat + bparsjj(ynj);
            pdfnorm2 = R::dnorm( -linpred2, 0.0, 1.0, 0L ); 
            cdfnorm2 = R::pnorm( -linpred2, 0.0, 1.0, 0L, 0L ); //  lower=false.
            Pi(0) = 1.0;
            Pi(1) = cdfnorm2;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 1) = pdfnorm2 * thetajj(temppar);
            //b-par
            dPidu(pj, 1) = pdfnorm2;
        }
        else if(ynj == mj_jj){
            linpred1 = sumat + bparsjj(ynj - 1);
            pdfnorm1 = R::dnorm( -linpred1, 0.0, 1.0, 0L ); 
            cdfnorm1 = R::pnorm( -linpred1, 0.0, 1.0, 0L, 0L );
            Pi(0) = cdfnorm1;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++) dPidu(temppar, 0) = pdfnorm1 * thetajj(temppar);
            //b-par
            dPidu(ynj - 2 + pj, 0) = pdfnorm1;
        }
        else{
            linpred1 = sumat + bparsjj(ynj - 1);
            linpred2 = sumat + bparsjj(ynj);
            pdfnorm1 = R::dnorm( -linpred1, 0.0, 1.0, 0L ); 
            cdfnorm1 = R::pnorm( -linpred1, 0.0, 1.0, 0L, 0L );
            pdfnorm2 = R::dnorm( -linpred2, 0.0, 1.0, 0L ); 
            cdfnorm2 = R::pnorm( -linpred2, 0.0, 1.0, 0L, 0L );
            Pi(0) = cdfnorm1;
            Pi(1) = cdfnorm2;
            //a-par
            for(arma::uword temppar = 0; temppar < pj; temppar++){
                dPidu(temppar, 0) = pdfnorm1 * thetajj(temppar);
                dPidu(temppar, 1) = pdfnorm2 * thetajj(temppar);
            }
            //b-par
            dPidu(ynj - 2 + pj, 0) = pdfnorm1;
            dPidu(ynj - 1 + pj, 1) = pdfnorm2;
        }
        
        //  SJ: Compute dPi / dtheta and d2Pi / d2theta
        for(arma::uword tt = 0; tt < 2; tt++){
            for(arma::uword uu = 0; uu < pj; uu++){
                dPidt(uu, tt) = (pdfnorm1 - pdfnorm2) / (Pi(0) - Pi(1)) * aparsjj(uu) ;
                for(arma::uword temppar = 0; temppar < nparj; temppar++){
                    d2Pidtdu(uu, temppar, tt) = aparsjj(uu) * (dPidu(temppar, tt) * (1.0 - Pi(tt)) - Pi(tt) * dPidu(temppar, tt));
                }
                d2Pidtdu(uu, uu, tt) += Pi(tt) * (1.0 - Pi(tt));
                for(arma::uword vv = 0; vv < pj; vv++){
                    d2Pidt2(uu, vv, tt) = aparsjj(uu) * dPidt(vv, tt) * (1.0 - 2.0 * Pi(tt));
                }
            }
        }
 
    }
    
    return(ret);
    
}

// [[Rcpp::export]]
Rcpp::List tabletolist(arma::mat estfixpars, arma::uword J, arma::vec mi, Rcpp::List model, std::vector<std::string> modeltype, arma::uword G){
    Rcpp::List apars(G);
    Rcpp::List bpars(G);
    arma::uword parindex = 0;
	arma::uvec npar(J);
	for(arma::uword g = 0; g < G; g++){
		Rcpp::List aparsg(J);
		Rcpp::List bparsg(J);
		apars(g) = aparsg;
		bpars(g) = bparsg;
	}
	for(arma::uword i = 0; i < J; i++){
		arma::uvec dimi = model(i);
		arma::uword pi = static_cast<uword>(dimi.n_elem);
		arma::uword npari = 0;
		arma::vec aparsi;
		arma::vec bparsi;
		for(arma::uword g = 0; g < G; g++){
			Rcpp::List aparsg = apars(g);
			Rcpp::List bparsg = bpars(g);
			if(modeltype[i] == "GPCM" || modeltype[i] == "GRM"){
				npari = pi + (mi(i) - 1);
				aparsi = zeros<vec>(pi);
				for(arma::uword j = 0; j < pi; j++){
					aparsi(j) = estfixpars(parindex + j, 1);
				}			
				//Add a zero to bpars for first parameter
				bparsi = zeros<vec>(mi(i));
				bparsi(0) = 0.0;
				bparsi(span(1, mi(i) - 1)) = estfixpars(span(parindex + pi, parindex + npari - 1), span(1));
				parindex += npari;
			}
			if(modeltype[i] == "NRM"){
				npari = pi * (mi(i) - 1) + (mi(i) - 1);
				aparsi = zeros<vec>(pi * (mi(i) - 1));
				for(arma::uword j = 0; j < pi * (mi(i) - 1); j++){
					aparsi(i) = estfixpars(parindex + j, 1);
				}				
				//Add a zero to bpars for first parameter
				bparsi = zeros<vec>(mi(i));
				bparsi(0) = 0.0;
				bparsi(span(1, mi(i) - 1)) = estfixpars(span(parindex + pi * (mi(i) - 1), parindex + npari - 1), span(1));
				parindex += npari;
			}
			if(modeltype[i] == "negbin" || modeltype[i] == "normal"){
				npari = pi + 2;
				aparsi = zeros<vec>(pi);
				for(arma::uword j = 0; j < pi; j++){
					aparsi(j) = estfixpars(parindex + j, 1);
				}				
				//Add a zero to bpars for first parameter
				bparsi = zeros<vec>(3);
				bparsi(0) = 0.0;
				bparsi(span(1, 2)) = estfixpars(span(parindex + pi, parindex + pi + 1), span(1));
				parindex += npari;
			}
			if(modeltype[i] == "poisson"){
				npari = pi + 1;
				aparsi = zeros<vec>(pi);
				for(arma::uword j = 0; j < pi; j++){
					aparsi(j) = estfixpars(parindex + j, 1);
				}				
				//Add a zero to bpars for first parameter
				bparsi = zeros<vec>(2);
				bparsi(0) = 0.0;
				bparsi(span(1, 1)) = estfixpars(span(parindex + pi, parindex + pi), span(1));
				parindex += npari;
			}
			aparsg(i) = aparsi;
			apars(g) = aparsg;
			bparsg(i) = bparsi;
			bpars(g) = bparsg;
		}
		npar(i) = npari;
	}
	return Rcpp::List::create(Rcpp::Named("apars") = apars,
                              Rcpp::Named("bpars") = bpars,
							  Rcpp::Named("npar") = npar);
}

// [[Rcpp::export]]
Rcpp::List mglogLGrad(arma::vec pars, arma::mat estfixpars, arma::mat y, arma::mat theta, arma::uword J, arma::vec mi, arma::uword p, Rcpp::List model, std::vector<std::string> modeltype, std::vector<std::string> link, arma::uword N, arma::mat covstruct, arma::uword G, arma::vec group, Rcpp::List filters, arma::mat X, std::string approx, arma::uword accuracy, arma::uvec npartype){
	//BA 2022-07-05: Change indices to j, k, l, m, n everywhere, and use z instead of t
    //BA 2022-07-02: New input argument arma::mat estfixpars;
    //BA 2022-07-02: Contains the non-zero parameters for each group and the indicator of the unique parameter each row corresponds to (NA for fixed parameters)
	//BA 2022-07-02: First column: Estimated (1) or not estimated (0)
	//BA 2022-07-02: Second column: Parameter value (estimate or constant)
	//BA 2022-07-02: Third column: Correspondence between unique parameters and estimated parameters
	//BA 2022-07-02: We append phi-pars to bpars for "negbin" and "normal"	
	//BA 2022-07-02: We need to know the total number of uniquely estimated parameters of each type in order to define the objects correctly.
	//BA 2022-07-02: Uniquepars object in R has information we need. Added a vector to this that indicates type of parameter: (item, regression, mean, variance, covariance)	
	//BA 2022-07-03: New input argument arma::uvec npartype
	//BA 2022-07-03: Entries indicate 1) the number of unique item parameters; 2) the number of unique regression parameters; 3) the number of unique mean parameters; 4) the number of unique variance parameters; 5) the number of unique covariance parameters
	
	//Define number of uniquely estimated parameters of each type
	//These determine the sizes of the gradient objects
	arma::uword nitpar = npartype(0);
	arma::uword nregpar = npartype(1);
	arma::uword nmupar = npartype(2);
	arma::uword nvarpar = npartype(3);
	arma::uword ncovpar = npartype(4);
	arma::uword nuniquepar = nitpar + nregpar + nmupar + nvarpar + ncovpar;

	//Initialize item parameter objects
	//We fill out these based on the input
    Rcpp::List apars(G);
    Rcpp::List bpars(G);
	
	//Define starting point for parameter index
    arma::uword parindex = 0;
	
	//This links the "active" parameters (estimated or fixed to constants) to the uniquely estimated parameters
	//For estimated parameters, vector indicates the position in the gradient vector (from 1 to nuniquepar). For fixed parameters, entry is 0. 
    arma::uvec fulltounique = conv_to<uvec>::from(estfixpars.col(2));

	//This prepares the item parameters into required lists per group, and specifies the number of uniquely estimated parameters per group
	Rcpp::List mylists = tabletolist(estfixpars, J, mi, model, modeltype, G);
	apars = mylists(0);
	bpars = mylists(1);
	arma::vec tmpobj = mylists(2);
	arma::uvec npar = conv_to<uvec>::from(tmpobj);

	//These are the numbers of "active" parameters of each type
	//BA 2022-07-05: (regression parameters not yet supported)	
	arma::uword nitpartot = G * sum(npar);
	arma::uword nregpartot = X.n_cols - 1;
	nregpartot = G * nregpartot * p;
	arma::uword nmupartot = G * p;
	arma::uword nvarpartot = G * p;
	arma::uword ncovpartot = G * p * (p - 1) / 2;
	
	//This is the total number of "active" parameters (the number of rows in estfixpars)
	arma::uword npartot = nitpartot + nregpartot + nmupartot + nvarpartot + ncovpartot;
	
	//Initialize means and covariance matrices
    arma::mat mu = zeros<mat>(p, G);
    arma::cube Sigma = zeros<cube>(p, p, G);
	
	//Rcout << "Run001";
	//BA 2022-08-25: We run through all the entries in the specific order of the input object
	//BA 2022-08-25: We update the index here to the point where the distribution parameters begin
	parindex = G * sum(npar);	
    for(arma::uword g = 0; g < G; g++){
		mu.col(g) = estfixpars(span(parindex, parindex + p - 1), span(1));
        parindex += p;			
    }
	//Rcout << "Run002";
	for(arma::uword g = 0; g < G; g++){
		Sigma.slice(g) = diagmat(estfixpars(span(parindex, parindex + p - 1), span(1)));
        parindex += p;
	}
	for(arma::uword g = 0; g < G; g++){
        //fill in covariance parameters
		for(arma::uword j = 0; j < (p - 1); j++){
			for(arma::uword k = j + 1; k < p; k++){
				Sigma(j, k, g) = estfixpars(parindex, 1);
				Sigma(k, j, g) = estfixpars(parindex, 1);
				parindex += 1;
			}
        }
    }
	
    //Define filters for 2nd-order Laplace
    //Input is the object 'filters', a list with entries:
    //uniqi3, uniqi4, Uniq3, Uniq4, UniqComb_3rd_4, UniqComb_3rd_6, UniqComb (0, 1, 2, 3, 4, 5, 6)
    //List w/ J entries, List w/ J entries, matrix, matrix, matrix, matrix, matrix 
	Rcpp::List uniqi3 = filters(0); // SJ: unique 3rd-order derivatives that we need to compute.
    Rcpp::List uniqi4 = filters(1); // SJ: unique 4th-order derivatives that we need to compute.
    arma::mat Uniq3 = filters(2);
    arma::mat Uniq4 = filters(3);
    arma::mat UniqComb_3rd_4 = filters(4);
    arma::mat UniqComb_3rd_6 = filters(5);
    arma::mat UniqComb = filters(6);
    
    arma::uword nUniqComb_3rd_4 = static_cast<uword>(UniqComb_3rd_4.n_rows);
    arma::uword nUniqComb_3rd_6 = static_cast<uword>(UniqComb_3rd_6.n_rows);
    arma::uword nUniqComb = static_cast<uword>(UniqComb.n_rows);
	
	//BA 2022-07-19: For some reason, the code below does not work.
	/*
    Rcpp::List uniqi3; // SJ: unique 3rd-order derivatives that we need to compute.
    Rcpp::List uniqi4; // SJ: unique 4th-order derivatives that we need to compute.
    arma::mat Uniq3;
    arma::mat Uniq4;
    arma::mat UniqComb_3rd_4;
    arma::mat UniqComb_3rd_6;
    arma::mat UniqComb;
    
    arma::uword nUniqComb_3rd_4;
    arma::uword nUniqComb_3rd_6;
    arma::uword nUniqComb;
	
	if(accuracy == 2){
		uniqi3 = filters(0);
		uniqi4 = filters(1);
		arma::mat Uniq3 = filters(2);
		arma::mat Uniq4 = filters(3);
		arma::mat UniqComb_3rd_4 = filters(4);
		arma::mat UniqComb_3rd_6 = filters(5);
		arma::mat UniqComb = filters(6);
		nUniqComb_3rd_4 = static_cast<uword>(UniqComb_3rd_4.n_rows);
		nUniqComb_3rd_6 = static_cast<uword>(UniqComb_3rd_6.n_rows);
		nUniqComb = static_cast<uword>(UniqComb.n_rows);
	}
	*/
    //Number of regression parameters
    arma::uword nbetapar = X.n_cols - 1;
    nbetapar = nbetapar * p;
		
    double pdouble = p;

    arma::mat betamat;	
    arma::mat mydhdb;
    arma::cube myd2hdtdb;
    if(G == 1){
        if(nregpartot > 0){
            betamat = zeros<mat>(p, nregpartot / p + 1);
            arma::uword betaparind = npartot - nregpartot;
            for(arma::uword tt = 0; tt < p; tt++){
                betamat(span(tt, tt), span(1, nregpartot / p)) = trans(pars(span(betaparind, betaparind + nregpartot / p - 1)));
                betaparind += nregpartot / p;
            }
            //Rcout << "Run001";
            for(arma::uword ss = 0; ss < p; ss++){
                arma::mat temp1 = X(span(0, N - 1), span(1, nregpartot / p));
                arma::vec temp2 = trans(betamat(span(ss, ss), span(1, nregpartot / p)));
                betamat(span(ss, ss), 0) = -mean(temp1 * temp2);
            }
            //Rcout << "Run002";
            mydhdb = dhdb(theta, Sigma.slice(0), betamat, X, p, N);
            //Rcout << "Run003";
            myd2hdtdb = d2hdtdb(theta, Sigma.slice(0), betamat, X, p, N);
            //Rcout << "Run004";
        }
    }
    
    //Rcout << "Run02";
    //Define inverses of covariance matrices
    arma::cube invSigma = zeros<cube>(p,p,G);
    for(arma::uword g = 0; g < G; g++) invSigma.slice(g) = inv(Sigma.slice(g));
    //Define constants
    arma::vec logsqrtdet2pisigma = zeros<vec>(G);
    for(arma::uword g = 0; g < G; g++){
        arma::mat tmpsig = 2.0 * M_PI * Sigma.slice(g);
        double templogsqrtdet2pisigma = det(tmpsig);
        templogsqrtdet2pisigma = sqrt(templogsqrtdet2pisigma);
        templogsqrtdet2pisigma = log(templogsqrtdet2pisigma);
        logsqrtdet2pisigma(g) = templogsqrtdet2pisigma;
    }
	
    //These are the output variables.
    //Output log-likelihood and gradient for each individual.
    arma::vec loglik = zeros<vec>(N);
    arma::mat gloglik = zeros<mat>(nuniquepar, N);
    
	//Define objects
    arma::vec addvec1;
    arma::vec addvec2;
	arma::uword muindex;
	arma::uword varindex;
	arma::uword covindex;
    for(arma::uword n = 0; n < N; n++){
		//BA 2022-07-05: Enable escaping from C++ to R prompt.
		Rcpp::checkUserInterrupt();
        //Prepare all the things we need.
        //For each item, we add to the entries of these objects.
        arma::uword ng = group(n);
        arma::vec mung = mu.col(ng);
        arma::mat invSigmang = invSigma.slice(ng);
        double logsqrtdet2pisigmag = logsqrtdet2pisigma(ng);
        double hn = 0.0;
        arma::vec dhn = zeros<vec>(p);
        arma::mat d2hn = zeros<mat>(p, p);
        arma::cube d3hn = zeros<cube>(p, p, p);
        arma::vec dhndu = zeros<vec>(nuniquepar);
        arma::mat d2hndu = zeros<mat>(p, nuniquepar);
        arma::cube d3hndu = zeros<cube>(p, p, nuniquepar);
        if(nregpartot > 0) mung = betamat * trans(X.row(n));
        arma::vec temptheta = trans(theta.row(n));
        arma::vec tempvec = temptheta - mung;
		//Rcout << "Run003";
		if(nmupar > 0){
            arma::vec dhdmu = zeros<vec>(p);
            arma::mat d2hdtdmu = zeros<mat>(p, p);
			for(arma::uword g = 0; g < G; g++){
				if(group(n) != g) continue;
				dhdmu(span(0, p - 1)) = -invSigma.slice(g) * tempvec;
				d2hdtdmu(span(0, p - 1), span(0, p - 1)) = -invSigma.slice(g);
			}
			//BA 2022-07-03: We loop over all mu-parameters in estfixpars, but only add estimated ones to the entries in dhndu and d2hndu
			for(arma::uword u = 0; u < p; u++){
				if(fulltounique(nitpartot + group(n) * p + u) == 0) continue;
				muindex = fulltounique(nitpartot + group(n) * p + u) - 1;
				dhndu(span(muindex, muindex)) += dhdmu(u);
				d2hndu(span(0, p - 1), span(muindex, muindex)) += d2hdtdmu.col(u);
				muindex += 1;
			}
        }
		
		//Rcout << "Run004";
        arma::vec dhds = zeros<vec>(p + p * (p - 1) / 2 );
        arma::mat d2hdtds = zeros<mat>(p, p + p * (p - 1) / 2);
        arma::cube d3hdt2ds = zeros<cube>(p, p, p + p * (p - 1) / 2);
        
        //Definition of checkmat1 and checkmat2 from covarianceprep()
        //arma::mat checkmat1 = covstruct;
		//BA 2022-07-04: We always go through all entries - will waste some time but makes code easier. covstruct input not needed, but can keep for now to avoid issues elsewhere.
        arma::mat checkmat2 = zeros<mat>(p, 2);
        checkmat2 = join_cols(checkmat2, covstruct);
        for(arma::uword j = 0; j < p; j++){
            checkmat2(j, 0) = j;
            checkmat2(j, 1) = j;
        }
        arma::mat checkmat;
        arma::uword s;
        arma::mat ident = eye(p, p);
        //Derivatives for covariance matrix parameters
		//Rcout << "Run005";
		for(arma::uword g = 0; g < G; g++){
            if(group(n) != g) continue;
            checkmat = checkmat2;
            arma::mat tempmat = -0.5 * (2.0 * invSigma.slice(g) - invSigma.slice(g) % ident -  2.0 * invSigma.slice(g) * tempvec * tempvec.t() * invSigma.slice(g) + invSigma.slice(g) * tempvec * tempvec.t() * invSigma.slice(g) % ident);
			arma::mat tempmat1MG;
            //Rcout << "Run101";			
			//Variance parameters
			s = 0;
			for(arma::uword u = 0; u < p; u++){
				arma::uword u1 = checkmat(u, 0);
                arma::uword u2 = checkmat(u, 1);
                tempmat1MG = zeros<mat>(p, p);
                dhds(s) = -tempmat(u1, u2);
                tempmat1MG(u1, u2) = 1.0;
                d2hdtds.col(s) = -invSigma.slice(g) * tempmat1MG * invSigma.slice(g) * tempvec;
                s += 1;
			}
			//Covariance parameters
			for(arma::uword u = p; u < p + p * (p - 1) / 2; u++){
				arma::uword u1 = checkmat(u, 0);
                arma::uword u2 = checkmat(u, 1);
				tempmat1MG = zeros<mat>(p, p);
				dhds(s) = -tempmat(u1, u2);
				tempmat1MG(u1, u2) = 1.0;
                tempmat1MG(u2, u1) = 1.0;
                d2hdtds.col(s) = -invSigma.slice(g) * tempmat1MG * invSigma.slice(g) * tempvec;
                s += 1;
            }
        }
	
		//New code
		//Rcout << "Run006";
		for(arma::uword g = 0; g < G; g++){
            if(group(n) != g) continue;
            checkmat = checkmat2;
            arma::vec tempvec = temptheta - mung;
			s = 0;
			for(arma::uword u = 0; u < p; u++){
				arma::uword u1 = checkmat(u, 0);
                arma::uword u2 = checkmat(u, 1);
                arma::mat tempmatMG2 = zeros<mat>(p, p);
				tempmatMG2(u1, u2) = 1.0;
                d3hdt2ds(span(0, p - 1), span(0, p - 1), span(s, s)) = - invSigma.slice(g) * tempmatMG2 * invSigma.slice(g);
                s += 1;
			}
			for(arma::uword u = p; u < p + p * (p - 1) / 2; u++){
				arma::uword u1 = checkmat(u, 0);
                arma::uword u2 = checkmat(u, 1);
                arma::mat tempmatMG2 = zeros<mat>(p, p);
				tempmatMG2(u1, u2) = 1.0;
                tempmatMG2(u2, u1) = 1.0;
                d3hdt2ds(span(0, p - 1), span(0, p - 1), span(s, s)) = - invSigma.slice(g) * tempmatMG2 * invSigma.slice(g);
                s += 1;
			}
		}
		
       // Rcout << "Run04";
		//BA 2022-07-04: Separate placement of variance/covariance derivatives since they are ordered separately by group
		//BA 2022-07-03: Link variance parameters between estfixpars and unique pars
		//BA 2022-07-03: Identify correct parameter from estfixpars, and link to unique parameter (if estimated)
		for(arma::uword u = 0; u < p; u++){
			if(fulltounique(nitpartot + nmupartot + group(n) * p + u) == 0) continue;
			varindex = fulltounique(nitpartot + nmupartot + group(n) * p + u) - 1;
			dhndu(span(varindex, varindex)) += dhds(u);
            d2hndu(span(0, p - 1), span(varindex, varindex)) += d2hdtds.col(u);
            d3hndu(span(0, p - 1), span(0, p - 1), span(varindex, varindex)) += d3hdt2ds.slice(u);
		}
		//if( n == 0 ){
		//	Rcout << "Index: " << varindex << std::endl;
        //}
		//BA 2022-07-03: Link covariance parameters between estfixpars and unique pars
		//BA 2022-07-03: Identify correct parameter from estfixpars, and link to unique parameter (if estimated)
		for(arma::uword u = p; u < p + p * (p - 1) / 2; u++){
			if(fulltounique(nitpartot + nmupartot + nvarpartot + group(n) * p * (p - 1) / 2 + u - p) == 0) continue;
			covindex = fulltounique(nitpartot + nmupartot + nvarpartot + group(n) * p * (p - 1) / 2 + u - p) - 1;
			dhndu(span(covindex, covindex)) += dhds(u);
            d2hndu(span(0, p - 1), span(covindex, covindex)) += d2hdtds.col(u);
            d3hndu(span(0, p - 1), span(0, p - 1), span(covindex, covindex)) += d3hdt2ds.slice(u);
		}
		//if( n == 0 ){
		//	Rcout << "Index: " << covindex << std::endl;
        //}
	
        //Derivatives for regression parameters (only for single group case right now)
		//BA 2022-07-04: This needs updating, broken order now
        if(nregpartot != 0){
            dhndu(span(nitpar, nitpar + nregpar - 1)) = trans(mydhdb.row(n));
            arma::mat myd2hdtdbn = myd2hdtdb(span(n, n), span(0, p - 1), span(0, nregpar - 1));
            d2hndu(span(0, p - 1), span(nitpar, nregpar - 1)) = -myd2hdtdbn;
        }
		
        double epsilon;
        arma::vec gepsilon = zeros<vec>(nuniquepar);
        arma::vec tepsilon = zeros<vec>(p);
        
        //Define object sizes for unique third- and fourth-order derivatives
        //For each item, we will add the unique entries to these objects
        //These will then be combined in the final calculation of the likelihood
        //The order is:  third- or fourth-order derivative, derivatives with respect to the latent variables, derivatives with respect to the unknown parameters
        arma::uword nuniq3all = static_cast<uword>(Uniq3.n_rows);
        arma::uword nuniq4all = static_cast<uword>(Uniq4.n_rows);
        arma::uword ncolsuniq34 = 1 + p + nuniquepar;
        arma::mat d3hdt3uniq = zeros<mat>(nuniq3all, ncolsuniq34);
        arma::mat d4hdt4uniq = zeros<mat>(nuniq4all, ncolsuniq34);
        
        //Multivariate normal distribution
        arma::mat tmpp = -0.5 * trans(temptheta - mung) * invSigmang * (temptheta - mung);
        hn += -(tmpp(0,0) - logsqrtdet2pisigmag);
        dhn += trans(invSigmang) * (temptheta - mung);
        d2hn += invSigmang;
		//BA 2022-07-02: Here, we need to keep track of the new order of things in the object estfixpars
		//BA 2022-07-03: We start with first item, which are the first parameters in estfixpars
        parindex = 0;
        arma::uword myjindex;
        Rcpp::List aparsg = apars(ng);
        Rcpp::List bparsg = bpars(ng);
        //--SJ: If we are using GRM with probit link, we need to compute a constant term for third order derivative. 
        double d3hconst = 0.0;
        
        // if ( n == 0 ){
        //     Rcout << "Cov dhn = " << dhn << std::endl;
        // }
        
        //Add unique stuff to objects that are later combined and weighted
        //Everything in this loop depends on the different IRT models
        for(arma::uword i = 0; i < J; i++){
            double midouble = mi(i);
			double ydouble = y(n, i);
			double sumPimi = 0.0;
			double phi = 1.0;
            arma::vec mival = regspace(1.0, midouble);
            //  SJ: p = dimension of latent variables, including those with zero slopes.
            //      pi = number of latent variables with non-zero slopes
            arma::uvec dimi = model(i);
            arma::uword pi = static_cast<uword>(dimi.n_elem);
			//BA 2022-07-02: Need to be item-specific
			//BA 2022-07-03: Updated, we have a new object "npari" which gives the number of *non-zero* parameters for each item (same number in each group)
			arma::uword npari = npar(i);
			parindex += group(n) * npari;
           // if( n == 0 ){
            //    Rcout << "item " << i << " pi=" << pi << " mi(i)=" << mi(i) << " npari=" << npari << std::endl;
			//	Rcout << "Index: " << parindex << std::endl;
            //}
            //Missing value handling, just skip to the next item. Missing coded as 9999.
			//BA 2022-07-02: Need to change below based on group, since we have ordered the item parameters by item by group
			//BA 2022-07-03: Should be fixed now
            if(y(n, i) == 9999){
                parindex += (G - group(n)) * npari;
                continue;
            }
            //dimi defines the non-zero 3rd order derivatives for item i
            arma::vec aparsi = aparsg(i);
            arma::vec bparsi = bparsg(i);
            arma::vec thetai = temptheta(dimi);
			double linpred = bparsi(1) + sum(aparsi % thetai);
			double toexp2 = 2.0 * linpred;
			double toexp3 = 3.0 * linpred;
			double explinpred = exp(linpred);
			double exp2linpred = exp(toexp2);
			double exp3linpred = exp(toexp3);
            arma::vec Pi;
            //  SJ: OBS! For GRM with probit link, we are computing logPi, not Pi.
			//	BA: Also, output is different between GRM and GPCM? Should be, for efficiency. Derivatives for GPCM and NRM have same structure, however.
            arma::mat dPidt;
            arma::cube d2Pidt2;
            arma::mat dPidu;
            arma::cube d2Pidtdu;
            // SJ: Added for new function (item-wise)
            arma::mat dlogPrdt;
            arma::mat d2logPrdt2;
            arma::cube d3logPrdt3;
            arma::mat dlogPrdu;
            arma::mat d2logPrdtdu;
            arma::cube d3logPrdt2du;
            //Add NRM to the functions below
            //  SJ: Compute the gradient of h with respect to item parameters
            arma::vec dhndui;
            //if( n == 0 ){
            //     Rcout << "  item " << i << " start. Item type is " << modeltype[i] << " link is " << link[i] << " thetai has dimension " << thetai.size() << std::endl;
            //}
			//BA: GPCM functions can be defined for (at least) NRM, too. 
			//BA 2022-06-29: dgjdu_string and d2gjdtdu_string need to handle fixed parameters (or we just need to ignore them somewhere)
			//BA 2022-07-03: We compute everything and ignore the zero ones later
            if(modeltype[i] == "GPCM"){
                Pi = gi(thetai, aparsi, bparsi, modeltype[i], mi(i), ydouble);
                dPidt = dgidz(thetai, aparsi, bparsi, modeltype[i], Pi, mi(i), pi, ydouble);
                d2Pidt2 = d2gidz2(thetai, aparsi, bparsi, modeltype[i], Pi, dPidt, mi(i), pi, ydouble);
                dPidu = dgidu(thetai, aparsi, bparsi, modeltype[i], Pi, mi(i), pi, ydouble);
                d2Pidtdu = d2gidzdu(thetai, aparsi, bparsi, modeltype[i], Pi, dPidt, dPidu, mi(i), pi, ydouble);

                double tempobj = Pi(y(n, i) - 1);
                sumPimi = sum(Pi % mival);
                hn += -log(tempobj);
                dhndui = -dPidu.col(y(n, i) - 1) / tempobj;
            }
            else if(modeltype[i] == "GRM"){
				//BA: This gives cumulative probs
                if(link[i] == "logit"){
                    //  SJ: Using arma::mat is expected to be faster. I am simply using Rcpp::List for simplicity.
					//	BA: Speed does seem to be negatively affected.
					// 	BA: What's the output for: 1) itemi_GRM, 2) itemi_GRM_quad, 3) itemi_GRM_3rd ?
					//	BA: New functions just convenience? I mean, still need item-if statement.
                    //Rcpp::List GRM_List = dgjdu_GRM_string(y(n, i), thetai, aparsi, bparsi, link[i], mi(i), pi, npari) ;
                    Rcpp::List GRM_List = item_GRM(y(n, i), thetai, aparsi, bparsi, link[i], mi(i), pi, npari) ;
                    //arma::vec Pi_temp = GRM_List["Pi"];
                    //Pi = Pi_temp;
                    //Pi = gj_GRM_logit(y(n, i), thetai, aparsi, bparsi, mi(i)) ;
					// 	BA: Does this assignment sequence do something special?
					//  BA: Goal is to remove the list and use "itemi_GRM_quad" only
                    arma::mat dPidt_temp = GRM_List["dPidt"];
                    dPidt = dPidt_temp; // SJ: I need dPidt.
                    arma::cube d2Pidt2_temp = GRM_List["d2Pidt2"];
                    d2Pidt2 = d2Pidt2_temp;
                    arma::mat dPidu_temp = GRM_List["dPidu"];
                    dPidu = dPidu_temp;
                    arma::cube d2Pidtdu_temp = GRM_List["d2Pidtdu"];
                    d2Pidtdu = d2Pidtdu_temp;

                    // SJ: I have used one function for 1st and 2nd order derivatives, and one function for 3rd order derivatives.
                    arma::mat GRMmat = item_GRM_quad(y(n, i), thetai, aparsi, bparsi, link[i], mi(i), pi, npari) ;
                    //if( n == 0 ){
                    //    Rcout << " itemi_GRM_quad OK";
                    //}
                    Pi = zeros<vec>(2);
                    Pi(0) = GRMmat(0,0);
                    Pi(1) = GRMmat(1,0);
                    hn += -log(Pi(0) - Pi(1));
                    dlogPrdt = GRMmat.submat(0, 1, pi - 1, 1);
                    dlogPrdu = GRMmat.submat(pi, 1, pi + npari -1, 1);
                    dhndui = -1.0 * dlogPrdu;
                    //if( n == 0 ){
                    //    Rcout << " Assign values OK";
                    //}
                    //ret.submat(0, 1, pi + npari -1, 1)
                    arma::cube GRM3rdcube = item_GRM_3rd(y(n, i), thetai, aparsi, bparsi, link[i], mi(i), pi, npari, Pi, true) ;
                    d3logPrdt3 = GRM3rdcube( span(0, pi - 1), span(0, pi - 1), span(0, pi - 1) );
                    d3logPrdt2du = GRM3rdcube( span(0, pi - 1), span(0, pi - 1), span(pi, pi + npari - 1) );
                    d2logPrdt2 = GRMmat.submat(0, 2, pi -1, 1 + pi);
                    d2logPrdtdu = GRMmat.submat(pi, 2, pi + npari -1, 1 + pi).t();

                }
                else if(link[i] == "probit"){
                    Pi = gj_GRM_probit(y(n, i), thetai, aparsi, bparsi, mi(i)) ;//  SJ: Not only Pr().
                    //   = Pivec( span(0,1) ) ;
                    dPidu = dgjdu_GRM_probit(y(n, i), thetai, Pi(4), Pi(5), mi(i), pi, npari);
                    dPidt = dgjdt_GRM_probit(aparsi, Pi(4), Pi(5), pi);
                    d2Pidt2 = d2gjd2t_GRM_probit(y(n, i), aparsi, Pi(2), Pi(3), Pi(4), Pi(5), mi(i), pi);
                   // if( n==0 ){
                    //    Rcout << "  ynj=" << y(n, i);
                    //}
                    d2Pidtdu = d2gjdtdu_GRM_probit(y(n, i), thetai, aparsi, Pi(2), Pi(3), Pi(4), Pi(5), mi(i), pi, npari);
                    //if( n==0 ){
                    //    Rcout << "  Finish d2logPidtdu" << std::endl;
                    //}
                    d3hconst = d3gjd3t_GRM_probit( y(n, i), aparsi, Pi(2), Pi(3), Pi(4), Pi(5), mi(i), pi );
                    
                    hn += -log(Pi(0) - Pi(1));
                    dhndui = -1.0 * (dPidu.col(0) - dPidu.col(1));
                }
            } else if(modeltype[i] == "negbin"){
				phi = bparsi(2);
				double negbintolog1 = 1.0 / phi + explinpred;
				double negbintolog2 = 1.0 + phi * explinpred;
				double todigam1 = ydouble + 1.0 / phi;
				double todigam2 = 1.0 / phi;
				double togammaf = ydouble + 1.0;
				//hn += -(ydouble * linpred - (ydouble + 1.0 / phi) * log(negbintolog) + ydouble * log(phi) + lgamma(todigam1) - lgamma(togammaf) - lgamma(todigam2));
				//hn += -(ydouble * linpred - (ydouble + 1.0 / phi) * log(negbintolog) + ydouble * log(phi) + lgamma(todigam1) - lgamma(togammaf) - lgamma(todigam2));
				/*
				hn += 
				-lgamma(ydouble + 1.0 / phi) 
				+ lgamma(ydouble + 1.0) 
				+ lgamma(ydouble + 1.0 / phi) 
				- ydouble * linpred 
				+ ydouble * log(1.0 / phi + explinpred) 
				+ (1.0 / phi) * log(1.0 + phi * explinpred);
				*/				
				hn += -lgamma(todigam1) + lgamma(togammaf) + lgamma(todigam2) - ydouble * linpred + ydouble * log(negbintolog1) + (1.0 / phi) * log(negbintolog2);
				//hn += -(ydouble * linpred - (ydouble + 1.0 / phi) * log(negbintolog) + ydouble * log(phi) + R::lgamma(todigam1) - R::lfactorial(togammaf) - R::lgamma(todigam2));
                dhndui = zeros<vec>(npari);
				for(arma::uword pari = 0; pari < pi; pari++){
					dhndui(pari) = -ydouble * thetai(pari) + (ydouble + 1.0 / phi) * (phi * explinpred) / (1.0 + phi * explinpred) * thetai(pari);
				}
				dhndui(pi) = -ydouble + (ydouble + 1.0 / phi) * (phi * explinpred) / (1.0 + phi * explinpred);
				//BA 2022-08-23: Need to check digamma function is actually available
				dhndui(pi + 1) = -log(negbintolog2) / (phi * phi) + (ydouble + 1.0 / phi) * (explinpred) / (1.0 + phi * explinpred) - ydouble / phi + R::digamma(todigam1) / (phi * phi) - R::digamma(todigam2) / (phi * phi);
				//if( n==0 ){
                //    Rcout << "dhndui done" << std::endl;
                //}
			} else if(modeltype[i] == "normal"){
				phi = bparsi(2);
				hn += -(ydouble * linpred - linpred * linpred / 2.0 ) / phi + ydouble * ydouble / (2.0 * phi) + (log2pi + log(phi)) / 2.0;
                dhndui = zeros<vec>(npari);
				for(arma::uword pari = 0; pari < pi; pari++){
					dhndui(pari) = thetai(pari) * (linpred - ydouble) / phi;
				}
				dhndui(pi) = (linpred - ydouble) / phi;
				//dhndui(pi + 1) = (phi - (linpred - ydouble) * (linpred - ydouble)) / (2.0 * phi * phi);
				dhndui(pi + 1) = ydouble * linpred / pow(phi, 2.0) - pow(linpred, 2.0) / (2.0 * pow(phi, 2.0)) - pow(ydouble, 2.0) / (2.0 * pow(phi, 2.0)) + 1.0 / (2.0 * phi);
			} else if(modeltype[i] == "poisson"){
				double yfac = rcpp_factorial(ydouble);
				hn += -log(pow(explinpred, ydouble) * exp(-explinpred) / yfac);
				//BA 2022-08-23: Need to implement
                dhndui = zeros<vec>(npari);
				for(arma::uword pari = 0; pari < pi; pari++){
					dhndui(pari) = thetai(pari) * (explinpred - ydouble);
				}
				dhndui(pi) = explinpred - ydouble;			
			}
			//BA 2022-07-04: Add NRM
			/*
			else if(modeltype[i] == "NRM"){
				//BA 2022-07-03: Need to implement
				hn += 0.0;
                dhndui = zeros<vec>(npari);
			} 
			*/
			//BA 2022-07-04: Can add gamma, beta, others here
			//else if(modeltype[i] == "gamma"){
			//	
			//}
			//else if(modeltype[i] == "beta"){
			//	
			//}
            //if( n==0 ){
            //    Rcout << "  item " << i << " ready.";
            //}
            
            //  SJ: Copy the gradient for item j into the gradient of all parameters
            //      This part does not depend on the item type.
            //apars
			//BA 2021-10-07: Seems like right place
			//BA 2022-06-29: Just copy derivatives for estimated pars?
			//BA 2022-07-02: Index over all parameters for the item, check if it is estimated or not, if it is not: skip to next parameter, if it is: find location and place right spot
			//And that's it.
            for(arma::uword pari = 0; pari < npari; pari++){
				//BA 2022-07-02: Need to account for some parameters not estimated
				if(fulltounique(parindex + pari) == 0) continue;
                myjindex = fulltounique(parindex + pari) - 1;
                dhndu(myjindex) += dhndui(pari);
            }
           // if( n==0 ){
            //    Rcout << "  item " << i << " assign ready.";
            //}
            //  SJ: The outer for-loop compute dh / dtheta (dhn)
            //  SJ: Why not compute the following loop in external modeltype specific functions?
            //      For the item parameters, we need dlogPr / dtheta
            //      What is addvec1??
			//	BA: addvec1 is just a temporary object that is added to the output objects
            //Rcout << "Run05";
            //Index variables
            //  SJ: If I understand the code correctly, index1 to index3 are used to find out which column of the complete apar we are talking about.
            //      Sometimes, we are only working with the nonzero entries in apar(i).
            //      But here, we need to work with the entire apar(i) vector.
            arma::uword index1 = 0;
            arma::uword index2;
            arma::uword index3;
            //This loop does repeated calculations for multidimensional models and can be improved. However, we need to fill out these objects in the end so changing it might not be worth it.
            for(auto j : dimi){
                index2 = 0;
                //posvec is defined as 1.0 for the value such that index1 is equal to the parameter value, i.e. the index1-th entry
				//BA 2022-06-29: npari needs to be the total number of non-zero pars, if we just copy the non-fixed pars later
                arma::vec posvec = zeros<vec>(npari);
                addvec1 = zeros<vec>(npari);
                //-- SJ: if(modeltype[i] == 2 || modeltype[i] == 4){
                if(modeltype[i] == "GPCM"){
                    dhn(j) += -aparsi(index1) * (ydouble - sumPimi);
                    posvec(index1) = 1.0;
                    for(arma::uword v = 0; v < npari; v++) addvec1(v) = aparsi(index1) * sum(trans(dPidu.row(v)) % mival);
                    addvec1 += -posvec * (ydouble - sumPimi);
                    //if( n == 0 ){
                    //    Rcout << "    item " << i << " dhn part = " << -aparsi(index1) * (ydouble - sumPimi) << std::endl;
                    //}
                } else if(modeltype[i] == "GRM"){
                    if(link[i] == "logit"){
                        dhn(j) += -dlogPrdt(index1,0);
                        //if( n == 0 ){
                        //    Rcout << "    item " << i << " dhn part = " << -dlogPrdt(index1,0) << std::endl;
                        //}
                        // SJ: Any better way to convert a row matrix to vector?
                        addvec1 = conv_to< arma::vec >::from( -1.0 * d2logPrdtdu.row(index1) ); 
                    }
                    else if(link[i] == "probit"){
                        dhn(j) += -aparsi(index1) * (Pi(4) - Pi(5));
                        for(arma::uword v = 0; v < npari; v++) addvec1(v) = aparsi(index1) * (dPidu(v, 0) + dPidu(v, 1));
                        addvec1(index1) += -(Pi(4) - Pi(5));
                    }
                } else if(modeltype[i] == "negbin"){
					dhn(j) += -ydouble * aparsi(index1) + (phi * ydouble + 1.0) * explinpred / (1.0 + phi * explinpred) * aparsi(index1);
					for(arma::uword v = 0; v < pi; v++){
						double myeq = (v == index1);
						addvec1(v) = -myeq * ydouble + (phi * ydouble + 1.0) * (explinpred) / (1.0 + phi * explinpred) * aparsi(index1) * thetai(v) - (phi * ydouble + 1.0) * (phi * exp2linpred) / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(index1) * thetai(v) + myeq * (phi * ydouble + 1.0) * explinpred / (1.0 + phi * explinpred) * (aparsi(index1) / aparsi(v));
					}
					addvec1(pi) = (phi * ydouble + 1.0) * explinpred / (1.0 + phi * explinpred) * aparsi(index1) - (phi * ydouble + 1.0) * (phi * exp2linpred) / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(index1);
					addvec1(pi + 1) = ydouble * (explinpred) / (1.0 + phi * explinpred) * aparsi(index1) - (phi * ydouble + 1.0) * exp2linpred / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(index1);
					//if(n == 0){
					//	Rcout << "d2hndu done" << std::endl;
					//}
				} else if(modeltype[i] == "normal"){
					dhn(j) += aparsi(index1) * (linpred - ydouble) / phi;
					for(arma::uword v = 0; v < pi; v++){
						double myeq = (v == index1);
						addvec1(v) = (aparsi(index1) * thetai(v) + (linpred - ydouble) * myeq) / phi;
					}
					addvec1(pi) = aparsi(index1) / phi;
					//addvec1(pi + 1) = -aparsi(index1) * (linpred - ydouble) / (phi * phi);
					addvec1(pi + 1) = ydouble * aparsi(index1) / pow(phi, 2.0) - linpred * aparsi(index1) / pow(phi, 2.0);
				} else if(modeltype[i] == "poisson"){
					dhn(j) += (explinpred - ydouble) * aparsi(index1);
					for(arma::uword v = 0; v < pi; v++){
						double myeq = (index1 == v);
						addvec1(v) = aparsi(index1) * explinpred * thetai(v)  + (explinpred - ydouble) * myeq;
					}
					addvec1(pi) = aparsi(index1) * explinpred;
				}
                //apars
				//BA 2022-07-03: Updated, we place the parameters in the correct place if estimated, and simply move to the next parameter if fixed to a constant
				for(arma::uword pari = 0; pari < npari; pari++){
					if(fulltounique(parindex + pari) == 0) continue;
					myjindex = fulltounique(parindex + pari) - 1;
					d2hndu(j, myjindex) += addvec1(pari);
				}
                //Rcout << "Run24";
                //if( n==0 ){
                //    Rcout << "  1st OK ";
               // }
                //  SJ: The inner for-loop computes d2h / d2theta (d2hn), the Hessian of h with respect to latent variables
                for(auto k : dimi){
                    index3 = 0;
                    addvec1 = zeros<vec>(npari);
                    //-- SJ: if(modeltype[i] == 2 || modeltype[i] == 4){
                    if(modeltype[i] == "GPCM"){
                        d2hn(j, k) += aparsi(index1) * sum(trans(dPidt.row(index2)) % mival);
                        arma::mat mytemp1 = d2Pidtdu(span(index2), span(0, npari - 1), span(0, mi(i) - 1));
                        for(arma::uword v = 0; v < npari; v++) addvec1(v) = aparsi(index1) * sum(trans(mytemp1.row(v)) % mival);
                        addvec1 += posvec * sum(trans(dPidt.row(index2)) % mival);
                    } else  if(modeltype[i] == "GRM"){
                        if(link[i] == "logit"){
                            d2hn(j, k) += -d2logPrdt2(index1, index2);
                            // SJ: I need to use arma::rowvec for subcube view. arma::vec does not work.
                            arma::rowvec direct_subset = d3logPrdt2du(arma::span(index1,index1), arma::span(index2,index2), arma::span::all);
                            addvec1 = conv_to< arma::vec >::from( -1.0 * direct_subset );
                        }
                        else if(link[i] == "probit"){
                            //  SJ:  Is this index correct???
							//	BA: Seems like this is incomplete? I haven't programmed the probit link so not sure what is going on here.
                            d2hn(j, k) += -d2Pidt2(index1,index2,0);
                            //  SJ:  How to make this addvec1 right??
                            //for(arma::uword v = 0; v < npari; v++) addvec1(v) = d2Pidtdu(index2, v, 0);
                            //addvec1(index1) += -d2Pidt2(i,j,0);
                        }
                    } else if(modeltype[i] == "negbin"){
						double topow  = 1.0 + phi * exp(linpred);
						d2hn(j, k) += (phi * ydouble + 1.0) * (explinpred) / pow(topow, 2.0) * aparsi(index1) * aparsi(index2);
						for(arma::uword v = 0; v < pi; v++){
							double myeq1 = (index1 == v);
							double myeq2 = (index2 == v);
							addvec1(v) = (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * aparsi(index1) * aparsi(index2) * thetai(v) - 2.0 * (phi * ydouble + 1.0) * (phi * exp2linpred) / pow(topow, 3.0) * aparsi(index1) * aparsi(index2) * thetai(v) + myeq1 * (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * aparsi(index2) + myeq2 * (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * aparsi(index1);
						}
						addvec1(pi) = (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * aparsi(index1) * aparsi(index2) - 2.0 * (phi * ydouble + 1.0) * (phi * exp2linpred) / pow(topow, 3.0) * aparsi(index1) * aparsi(index2);
						addvec1(pi + 1) = ydouble * explinpred / pow(topow, 2.0) * aparsi(index1) * aparsi(index2) - 2.0 * (phi * ydouble + 1.0) * exp2linpred / pow(topow, 3.0) * aparsi(index1) * aparsi(index2);
						//if(n == 0){
						//	Rcout << "d3hndu done" << std::endl;
						//}
					} else if(modeltype[i] == "normal"){
						d2hn(j, k) += aparsi(index1) * aparsi(index2) / phi;
						for(arma::uword v = 0; v < pi; v++){
							double myeq1 = (index1 == v);
							double myeq2 = (index2 == v);
							addvec1(v) = (aparsi(index2) * myeq1 + aparsi(index1) * myeq2) / phi;
						}
						//Derivative with respect to the intercept is zero
						addvec1(pi + 1) = -aparsi(index1) * aparsi(index2) / pow(phi, 2.0);
					} else if(modeltype[i] == "poisson"){
						d2hn(j, k) += explinpred * aparsi(index1) * aparsi(index2);
						for(arma::uword v = 0; v < pi; v++){
							double myeq1 = (index1 == v);
							double myeq2 = (index2 == v);
							addvec1(v) = aparsi(index1) * aparsi(index2) * explinpred * thetai(v) + explinpred * (aparsi(index2) * myeq1 + aparsi(index1) * myeq2);
						}
						addvec1(pi) = aparsi(index1) * aparsi(index2) * explinpred;
					}
					//BA 2022-07-04: Add NRM,
                    //if( n==0 ){
                    //    Rcout << "  2nd Assign ";
                    //}
                    //apars
					//BA 2021-10-07: Seems like right place
					//BA 2022-06-29: Just copy derivatives for estimated pars?
					//BA 2022-07-03: Updated
					for(arma::uword pari = 0; pari < npari; pari++){
						if(fulltounique(parindex + pari) == 0) continue;
						myjindex = fulltounique(parindex + pari) - 1;
						d3hndu(j, k, myjindex) += addvec1(pari);
					}
                    //if( n==0 ){
                    //    Rcout << "  2nd OK ";
                    //}
                    // SJ: the last for-loop to compute d3h / d3theta (d3hn)
					// Third order and higher derivatives are zero for modeltype[i] == "normal"
                    for(auto l : dimi){
                        //-- SJ: if(modeltype[i] == 2 || modeltype[i] == 4){
                        if(modeltype[i] == "GPCM"){
                            arma::vec mytemp2 = d2Pidt2(span(index2), span(index3), span(0, mi(i) - 1));
                            d3hn(j, k, l) += aparsi(index1) * sum(mytemp2 % mival);
                        } else if(modeltype[i] == "GRM"){
                            if(link[i] == "logit"){
                                d3hn(j, k, l) += -d3logPrdt3(index1, index2, index3);
                            }
                            else if(link[i] == "probit"){
                                d3hn(j, k, l) += d3hconst * aparsi(index1) * aparsi(index2) * aparsi(index3) ;
                            }
                        } else if(modeltype[i] == "negbin"){
							d3hn(j, k, l) += -(phi * ydouble + 1.0) * explinpred * (phi * explinpred - 1.0) / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(index1) * aparsi(index2) * aparsi(index3); 
						} else if(modeltype[i] == "poisson"){
							d3hn(j, k, l) += explinpred * aparsi(index1) * aparsi(index2) * aparsi(index3);
						}
						//BA 2022-07-04: Add NRM
                        //if( n==0 ){
                        //    Rcout << "  3rd OK ";
                        //}
                        index3 += 1;
                    }
                    index2 += 1;
                }
                index1 += 1;
            }
            //if( n==0 ){
            //   Rcout << " index1=" << index1 << " Further loops ready.";
           // }
            
            //  If we are working with second-order Laplace approximation
            if( approx == "Laplace" && accuracy == 2 ){
                //Rcout << "Run06";
                double d3hidt3uniq = 0.0;
                double d4hidt4uniq = 0.0;
                // SJ: When we have NRM with probit link, some higher-order derivatives are simply constants * parameter-specific terms. 
                //     Hence, we can compute constant terms only once.
                double d4hconst = 0.0;
                double d5hconst = 0.0;
                if(modeltype[i] == "GRM" && link[i] == "probit"){
                    d4hconst = d4gjd4t_GRM_probit(y(n, i), aparsi, Pi(2), Pi(3), Pi(4), Pi(5), mi(i), pi);
                    d5hconst = d5gjd5t_GRM_probit(y(n, i), aparsi, Pi(2), Pi(3), Pi(4), Pi(5), mi(i), pi);
                }
                // if( n==0 ){
                //     Rcout << " Constants ready." << std::endl;
                // }
                // SJ: uniqi3h contains the entries of 3rd-order derivatives needed for second-order Laplace loops. 
                arma::uvec uniqi3h = uniqi3(i);
                //if( n==0 ){
                //     Rcout << " nuniq3all=" << nuniq3all << " ncolsuniq34=" << ncolsuniq34 << std::endl;
                // }
                
                //  SJ: What we have here is functions used in for(auto h : uniqi3h).
                //  An alternative way is to put for(auto h : uniqi3h) in a function, stil one func for one item.
                //  The code in the main function is neater, but tempj, tempk, templ part will be repated in each function.
                //  A stand-alone question: higher derivatievs of the latent regression part? 
				//	BA: higher-order derivatives are zero with the model I programmed.
                for(auto h : uniqi3h){
                    // SJ: The only object needs to be returned is arma::mat d3hdt3uniq = zeros<mat>(nuniq3all, ncolsuniq34);
                    //Which dimensions do we need?
                    // SJ: d3h / dtheta(j) dtheta(k) dtheta(l). (j,k,l) is the index of the full dimension.
                    arma::uword j = Uniq3(h, 0);
                    arma::uword k = Uniq3(h, 1);
                    arma::uword l = Uniq3(h, 2);
                    //if( n == 0){
                    //    Rcout << " h=" << h << std::endl;
                    //}
                    
                    //Go from full dimensions (1, ..., p) to the dimensions for the item (1, ..., pit)
                    // SJ: (j, k, l) is the index of the full dimension.
                    //     (tempj, tempk, templ) is the index of the reduced dimension, i.e. only those with nonzero apars.
                    arma::uword tempj = 99;
                    arma::uword tempk = 99;
                    arma::uword templ = 99;
                    for(arma::uword ii = 0; ii < pi; ii++){
                        if(dimi(ii) == j) tempj = ii;
                        if(dimi(ii) == k) tempk = ii;
                        if(dimi(ii) == l) templ = ii;
                    }
                    //Rcout << "Run06A";
                    // SJ: Under which conditions will temp index not found? 
                    if(tempj == 99 || tempk == 99 || templ == 99) continue;
                    //if(tempj == 99) continue;
                    //if(tempk == 99) continue;
                    //if(templ == 99) continue;
                    
                    arma::vec sumeach;
                    arma::vec mival;
                    //Write functions to calculate d3Pidt2kldu and d3Pidt3klm?
                    //Input: d3Pidt2kldu - Pi, dPidt, d2Pidt2, dPidu, d2Pidtdu, mival, aparsi, const1, tempk, templ
                    arma::mat d3Pidt2kldu = zeros<mat>(npari, mi(i)); // SJ: Do we need this in GRM with probit link? BA: You tell me! :)
                    arma::vec d4hidt3duuniq = zeros<vec>(npari); // SJ: Is this d4hi / dtheta^3 dpar ? BA: Yes.
                    arma::vec d4hidt4jklm = zeros<vec>(pi);
                    if(modeltype[i] == "GPCM"){
                        sumeach = zeros<vec>(mi(i));
                        mival = regspace(1.0, midouble);
                        for(arma::uword u = 0; u < npari; u++){
                            double myind = tempk == u;
                            arma::vec tempvec = d2Pidtdu.tube(templ,u);
                            for(arma::uword cat = 0; cat < mi(i); cat++){	
                                d3Pidt2kldu(u, cat) = myind * d2Pidt2(tempk, templ, cat) / aparsi(tempk) + d2Pidtdu(templ, u, cat) * mival(cat) * aparsi(tempk) - d2Pidtdu(templ, u, cat) * aparsi(tempk) * sumPimi - dPidt(templ, cat) * aparsi(tempk) * sum(trans(dPidu.row(u)) % mival) - dPidu(u, cat) * aparsi(tempk) * sum(trans(dPidt.row(templ)) % mival) - Pi(cat) * aparsi(tempk) * sum(tempvec % mival);
                            }
                        }
                        arma::vec tempvec = d2Pidt2.tube(tempk,templ);
                        d3hidt3uniq = aparsi(tempj) * sum(tempvec % mival);
                        for(arma::uword u = 0; u < npari; u++){
                            double myind = tempj == u;
                            d4hidt3duuniq(u) = myind * sum(tempvec % mival) + aparsi(tempj) * sum(trans(d3Pidt2kldu.row(u)) % mival);
                        }				
                        for(arma::uword m = 0; m < pi; m++){
                            arma::vec d3Pidt3klm = zeros<vec>(mi(i));
                            arma::vec d2Pidt2lm = d2Pidt2.tube(templ, m);
                            //Define constants here
                            for(arma::uword cat = 0; cat < mi(i); cat++){
                                d3Pidt3klm(cat) = d2Pidt2(templ, m, cat) * mival(cat) * aparsi(tempk) - d2Pidt2(templ, m, cat) * aparsi(tempk) * sumPimi - dPidt(templ, cat) * aparsi(tempk) * sum(trans(dPidt.row(m)) % mival) - dPidt(m, cat) * aparsi(tempk) * sum(trans(dPidt.row(templ)) % mival) - Pi(cat) * aparsi(tempk) * sum(d2Pidt2lm % mival);
                            }
                            d4hidt4jklm(m) = aparsi(tempj) * sum(d3Pidt3klm % mival);
                        }
                    } else if(modeltype[i] == "GRM"){
                        if(link[i] == "logit"){
                            // SJ: We only need to return d3Pidt2kldu, d4hidt3duuniq, d4hidt4jklm, d3hidt3uniq
                            for(arma::uword u = 0; u < npari; u++){
                                d3Pidt2kldu(u, 0) = aparsi(tempk) * (d2Pidtdu(templ, u, 0) * (1.0 - 2.0 * Pi(0)) - 2.0 * dPidt(templ, 0) * dPidu(u, 0));
                                d3Pidt2kldu(u, 1) = aparsi(tempk) * (d2Pidtdu(templ, u, 1) * (1.0 - 2.0 * Pi(1)) - 2.0 * dPidt(templ, 1) * dPidu(u, 1));
                            }
                            d3Pidt2kldu(tempk, 0) += dPidt(templ, 0) * (1.0 - 2.0 * Pi(0));
                            d3Pidt2kldu(tempk, 1) += dPidt(templ, 1) * (1.0 - 2.0 * Pi(1));
                            for(arma::uword u = 0; u < npari; u++) d4hidt3duuniq(u) = aparsi(tempj) * (d3Pidt2kldu(u, 0) + d3Pidt2kldu(u, 1));
                            d4hidt3duuniq(tempj) += (d2Pidt2(tempk, templ, 0) + d2Pidt2(tempk, templ, 1));
                            for(arma::uword m = 0; m < pi; m++){					
                                arma::vec d3Pidt3klm = zeros<vec>(2);
                                d3Pidt3klm(0) = aparsi(tempk) * (d2Pidt2(templ, m, 0) * (1.0 - 2.0 * Pi(0)) - 2.0 * dPidt(templ, 0) * dPidt(m, 0));
                                d3Pidt3klm(1) = aparsi(tempk) * (d2Pidt2(templ, m, 1) * (1.0 - 2.0 * Pi(1)) - 2.0 * dPidt(templ, 1) * dPidt(m, 1));
                                d4hidt4jklm(m) = aparsi(tempj) * (d3Pidt3klm(0) + d3Pidt3klm(1));
                            }
                            d3hidt3uniq = aparsi(tempj) * (d2Pidt2(tempk, templ, 0) + d2Pidt2(tempk, templ, 1));
                        }
                        else if(link[i] == "probit"){
                            // for(arma::uword u = 0; u < npari; u++){
                            //     d3Pidt2kldu(u, 0) = aparsi(tempk) * (d2Pidtdu(templ, u, 0) * (1.0 - 2.0 * Pi(0)) - 2.0 * dPidt(templ, 0) * dPidu(u, 0));
                            //     d3Pidt2kldu(u, 1) = aparsi(tempk) * (d2Pidtdu(templ, u, 1) * (1.0 - 2.0 * Pi(1)) - 2.0 * dPidt(templ, 1) * dPidu(u, 1));
                            // }
                            // d3Pidt2kldu(tempk, 0) += dPidt(templ, 0) * (1.0 - 2.0 * Pi(0));
                            // d3Pidt2kldu(tempk, 1) += dPidt(templ, 1) * (1.0 - 2.0 * Pi(1));
                            // for(arma::uword u = 0; u < npari; u++) d4hidt3duuniq(u) = aparsi(tempj) * (d3Pidt2kldu(u, 0) + d3Pidt2kldu(u, 1));
                            // d4hidt3duuniq(tempj) += (d2Pidt2(tempk, templ, 0) + d2Pidt2(tempk, templ, 1));
                            // SJ: apar
                            for(arma::uword u = 0; u < pi; u++){
                                d4hidt3duuniq(u) = d4hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * thetai(u) + (aparsi(tempk) * aparsi(templ) + aparsi(tempj) * aparsi(templ) + aparsi(tempj) * aparsi(tempk)) * d3hconst; // SJ: We need thetai to be the reduced dimension theta.
                            }
                            // SJ: bpar
                            if(y(n, j) == 1){
                                d4hidt3duuniq(pi) += d4hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ;
                            }
                            else if(y(n, j) == mi(i)){
                                d4hidt3duuniq(npari-1) += d4hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ; 
                            }
                            else{
                                d4hidt3duuniq(pi + y(n, j) - 2) += d4hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ;
                                d4hidt3duuniq(pi + y(n, j) - 1) += d4hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ;
                            }
                            
                            d4hidt4jklm = d4hconst * aparsi;
                            d3hidt3uniq = d3hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ;
                        }
                    } else if(modeltype[i] == "negbin"){
						d3hidt3uniq = -(phi * ydouble + 1.0) * explinpred * (phi * explinpred - 1.0) / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(tempj) * aparsi(tempk) * aparsi(templ);
						double topow = 1.0 + phi * explinpred;
						for(arma::uword u = 0; u < pi; u++){
							double myeq1 = (tempj == u);
							double myeq2 = (tempk == u);
							double myeq3 = (templ == u);
                            d4hidt3duuniq(u) = -(phi * ydouble + 1.0) * (explinpred * (phi * explinpred - 1.0)) / pow(topow, 3.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * thetai(u) - (phi * ydouble + 1.0) * (phi * exp2linpred) / pow(topow, 3.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * thetai(u) + 3.0 * (phi * ydouble + 1.0) * (phi * exp2linpred * (phi * explinpred - 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * thetai(u) - myeq1 * (phi * ydouble + 1.0) * (explinpred * (phi * explinpred - 1.0)) / pow(topow, 3.0) * aparsi(tempk) * aparsi(templ) - myeq2 * (phi * ydouble + 1.0) * (explinpred * (phi * explinpred - 1.0)) / pow(topow, 3.0) * aparsi(tempj) * aparsi(templ) - myeq3 * (phi * ydouble + 1.0) * (explinpred * (phi * explinpred - 1.0)) / pow(topow, 3.0) * aparsi(tempj) * aparsi(tempk);
                        }
						d4hidt3duuniq(pi) = - (phi * ydouble + 1.0) * (explinpred * (phi * explinpred - 1.0)) / pow(topow, 3.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) - (phi * ydouble + 1.0) * (phi * exp2linpred) / pow(topow, 3.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) + 3.0 * (phi * ydouble + 1.0) * (phi * exp2linpred * (phi * explinpred - 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ);						
						d4hidt3duuniq(pi + 1) = - ydouble * (explinpred * (phi * explinpred - 1.0)) / pow(topow, 3.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) - (phi * ydouble + 1.0) * (exp2linpred) / pow(topow, 3.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) + 3.0 * (phi * ydouble + 1.0) * (exp2linpred * (phi * explinpred - 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ);
						for(arma::uword m = 0; m < pi; m++){
                            d4hidt4jklm(m) = (phi * ydouble + 1.0) * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred) * (1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(m);
                        }
						//if(n == 0){
						//	Rcout << "d4hndu done" << std::endl;
						//}
					} else if(modeltype[i] == "poisson"){
                        d3hidt3uniq = explinpred * aparsi(tempj) * aparsi(tempk) * aparsi(templ);
                        for(arma::uword u = 0; u < pi; u++){
							double myeq1 = (tempj == u);
							double myeq2 = (tempk == u);
							double myeq3 = (templ == u);
                            d4hidt3duuniq(u) = aparsi(tempj) * aparsi(tempk) * aparsi(templ) * explinpred * thetai(u) + explinpred * (aparsi(tempk) * aparsi(templ) * myeq1 + aparsi(tempj) * aparsi(templ) * myeq2  + aparsi(tempj) * aparsi(tempk) * myeq3);
                        }
						d4hidt3duuniq(pi) = explinpred * aparsi(tempj) * aparsi(tempk) * aparsi(templ);
                        for(arma::uword m = 0; m < pi; m++){
                            d4hidt4jklm(m) = explinpred * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(m);
                        }
					}
					//BA 2022-08-23: Add NRM, negbin
                    //Rcout << "Run06C";
                    //Unique third order derivatives
                    //  SJ: Instead, we can use one function for each h to return the h-th row of d3hdt3uniq
					//	BA: If making a function doesn't add much overhead it's probably worth doing to de-clutter the code if nothing else.
					// 	BA: Don't follow this if statement. (This messes up GPCM? Commented out if/else.)
                    //if( i != 9 || modeltype[i] == "GPCM"){
                        d3hdt3uniq(h, 0) += d3hidt3uniq;
                        //Fourth derivatives with respect to latent variables
                        for(arma::uword ii = 0; ii < pi; ii++){
                            arma::uword myentry = dimi(ii);
                            d3hdt3uniq(h, 1 + myentry) += d4hidt4jklm(ii);
                        }			
                        //Rcout << "Run06D";				
                        //Derivatives with respect to item parameters
                        //apars
						//BA 2022-06-29: Just copy derivatives for estimated pars?
						//BA 2022-07-03: Updated
						for(arma::uword pari = 0; pari < npari; pari++){
							if(fulltounique(parindex + pari) == 0) continue;
							myjindex = fulltounique(parindex + pari) - 1;
							d3hdt3uniq(h, 1 + p + myjindex) += d4hidt3duuniq(pari);
						}
                    //} else {
                    //    arma::mat jsb = itemi_GRM_lap2_V3( y(n, i), thetai, aparsi, link[i], mi(i), pi, npari, p, ncolsuniq34, tempj, tempk, templ, Pi, dimi, parindex, fulltounique ) ;
                    //    d3hdt3uniq.row(h) += jsb;
                        //if( n == 0 ){
                        //    Rcout << " d3hdt3uniq(h,0) = " << d3hidt3uniq << " jsb(0,0)=" << jsb(0,0) << std::endl;
                        //    Rcout << " d4hidt4jklm = " << d4hidt4jklm.t() << std::endl;
                        //    Rcout << " jsb=" << jsb.submat(0,1,0,pi) << std::endl;
                        //    Rcout << " d4hidt3duuniq = " << d4hidt3duuniq.t() << std::endl;
                        //    Rcout << " jsb=" << jsb.submat(0,pi+1,0,ncolsuniq34-1) << std::endl;
                        //}
                    //}
                    
                    //d3hdt3uniq.row(h) += 
                   
                    //Rcout << "Run06F";				
                }
                 //if( n==0 ){
                //     Rcout << " 3rd-order ready." << std::endl;
                // }

                //Rcout << "Run07";
                //  SJ: We can do a similar thing to d4hdt4uniq.
                //      Either one function for each h or for(h) in the function.
				//BA 2022-07-04: Make h4-function, which has inputs: uniqi4h, dimi, pi, p, npari, mi, Pi, dPidt, d2Pidt2, sumPimi, aparsi, mival, dPidu, d2Pidtdu
				//Output: d4hidt4uniq, d5hidt5uniq, d5hidt4duuniq
                arma::uvec uniqi4h = uniqi4(i);
                for(auto h : uniqi4h){
                    //Which dimensions do we need?
                    arma::uword j = Uniq4(h, 0);
                    arma::uword k = Uniq4(h, 1);
                    arma::uword l = Uniq4(h, 2);
                    arma::uword m = Uniq4(h, 3);
                    
                    //Go from full dimensions (1, ..., p) to the dimensions for the item (1, ..., pit)
                    arma::uword tempj = 99;
                    arma::uword tempk = 99;
                    arma::uword templ = 99;
                    arma::uword tempm = 99;
                    for(arma::uword ii = 0; ii < pi; ii++){
                        if(dimi(ii) == j) tempj = ii;
                        if(dimi(ii) == k) tempk = ii;
                        if(dimi(ii) == l) templ = ii;
                        if(dimi(ii) == m) tempm = ii;
                    }
                    if(tempj == 99) continue;
                    if(tempk == 99) continue;
                    if(templ == 99) continue;
                    if(tempm == 99) continue;
                    arma::mat d3Pidt2lmdu;
                    arma::mat d4Pidt3klmdu;
                    arma::vec d5hidt4duuniq;
                    arma::mat d4Pidt4klmn;
                    arma::vec d5hidt5uniq = zeros<vec>(p);
                    
					//	BA: This code is pretty crazy.. Probably need to check that this works with all types of model structures
                    if(modeltype[i] == "GPCM"){
                        d3Pidt2lmdu = zeros<mat>(npari, mi(i));
                        d4Pidt3klmdu = zeros<mat>(npari, mi(i));
                        d5hidt4duuniq = zeros<vec>(npari);
                        d4Pidt4klmn = zeros<mat>(p, mi(i));
                        arma::vec d3Pidt3klm = zeros<vec>(mi(i));
                        arma::vec d2Pidt2lm = d2Pidt2.tube(templ, tempm);
                        //Define constants here
                        for(arma::uword cat = 0; cat < mi(i); cat++) d3Pidt3klm(cat) = d2Pidt2(templ, tempm, cat) * aparsi(tempk) * mival(cat) - d2Pidt2(templ, tempm, cat) * aparsi(tempk) * sumPimi - dPidt(templ, cat) * aparsi(tempk) * sum(trans(dPidt.row(tempm)) % mival) - dPidt(tempm, cat) * aparsi(tempk) * sum(trans(dPidt.row(templ)) % mival) - Pi(cat) * aparsi(tempk) * sum(d2Pidt2lm % mival);
                        //Rcout << "Run07A";
                        d4hidt4uniq = aparsi(tempj) * sum(d3Pidt3klm % mival);
                        for(arma::uword u = 0; u < npari; u++){
                            double myind = templ == u;
                            arma::vec d2Pidtdum = d2Pidtdu.tube(tempm, u);
                            //Define constants here
                            for(arma::uword cat = 0; cat < mi(i); cat++) d3Pidt2lmdu(u, cat) = myind * d2Pidt2(templ, tempm, cat) / aparsi(templ) + aparsi(templ) * d2Pidtdu(tempm, u, cat) * mival(cat)	- aparsi(templ) * d2Pidtdu(tempm, u, cat) * sumPimi - aparsi(templ) * dPidt(tempm, cat) * sum(trans(dPidu.row(u)) % mival) - aparsi(templ) * dPidu(u, cat) * sum(trans(dPidt.row(tempm)) % mival) - aparsi(templ) * Pi(cat) * sum(d2Pidtdum % mival);	
                        }
                        //Rcout << "Run07B";
                        for(arma::uword u = 0; u < npari; u++){
                            double myind = tempk == u;
                            arma::vec tempvec1 = d2Pidtdu.tube(tempm, u);
                            arma::vec tempvec2 = d2Pidtdu.tube(templ, u);
                            arma::vec tempvec3 = d2Pidt2.tube(templ, tempm);
                            for(arma::uword cat = 0; cat < mi(i); cat++) d4Pidt3klmdu(u, cat) = myind * d3Pidt3klm(cat) / aparsi(tempk) + aparsi(tempk) * d3Pidt2lmdu(u, cat) * mival(cat) - aparsi(tempk) * d3Pidt2lmdu(u, cat) * sumPimi - aparsi(tempk) * d2Pidt2(templ, tempm, cat) * sum(trans(dPidu.row(u)) % mival) - aparsi(tempk) * d2Pidtdu(templ, u, cat) * sum(trans(dPidt.row(tempm)) % mival) - aparsi(tempk) * dPidt(templ, cat) * sum(tempvec1 % mival) - aparsi(tempk) * d2Pidtdu(tempm, u, cat) * sum(trans(dPidt.row(templ)) % mival) - aparsi(tempk) * dPidt(tempm, cat) * sum(tempvec2 % mival) - aparsi(tempk) * dPidu(u, cat) * sum(tempvec3 % mival) - aparsi(tempk) * Pi(cat) * sum(trans(d3Pidt2lmdu.row(u)) % mival);
                        }
                        //Rcout << "Run07C";
                        for(arma::uword u = 0; u < npari; u++){
                            double myind = tempj == u;
                            d5hidt4duuniq(u) = myind * sum(d3Pidt3klm % mival) + aparsi(tempj) * sum(trans(d4Pidt3klmdu.row(u)) % mival);
                        }
                        for(arma::uword kk = 0; kk < p; kk++){
                            arma::uword tempn = 99;
                            for(arma::uword tt = 0; tt < pi; tt++){
                                if(dimi(tt) == kk) tempn = tt; 
                            }
                            if(tempn == 99) continue;
                            arma::vec d3Pidt3lmn = zeros<vec>(mi(i));
                            arma::vec tempvec1 = d2Pidt2.tube(tempm,tempn);
                            arma::vec tempvec2 = d2Pidt2.tube(templ,tempn);
                            arma::vec tempvec3 = d2Pidt2.tube(templ,tempm);
                            //Cross-check with d3gjdt3 (checks out)
                            for(arma::uword cat = 0; cat < mi(i); cat++) d3Pidt3lmn(cat) = aparsi(templ) * d2Pidt2(tempm, tempn, cat) * mival(cat) - aparsi(templ) * d2Pidt2(tempm, tempn, cat) * sumPimi - aparsi(templ) * dPidt(tempm, cat) * sum(trans(dPidt.row(tempn)) % mival) - aparsi(templ) * dPidt(tempn, cat) * sum(trans(dPidt.row(tempm)) % mival) - aparsi(templ) * Pi(cat) * sum(tempvec1 % mival);
                            //Cross-check with d4gjdt4 (checks out)
                            for(arma::uword cat = 0; cat < mi(i); cat++) d4Pidt4klmn(kk, cat) = aparsi(tempk) * d3Pidt3lmn(cat) * mival(cat) - aparsi(tempk) * d3Pidt3lmn(cat) * sumPimi - aparsi(tempk) * d2Pidt2(templ, tempm, cat) * sum(trans(dPidt.row(tempn)) % mival) - aparsi(tempk) * d2Pidt2(templ, tempn, cat) * sum(trans(dPidt.row(tempm)) % mival) - aparsi(tempk) * dPidt(templ, cat) * sum(tempvec1 % mival) - aparsi(tempk) * d2Pidt2(tempm, tempn, cat) * sum(trans(dPidt.row(templ)) % mival) - aparsi(tempk) * dPidt(tempm, cat) * sum(tempvec2 % mival) - aparsi(tempk) * dPidt(tempn, cat) * sum(tempvec3 % mival) - aparsi(tempk) * Pi(cat) * sum(d3Pidt3lmn % mival);
                            d5hidt5uniq(kk) += aparsi(tempj) * sum(trans(d4Pidt4klmn.row(kk)) % mival);
                        }
                        //Rcout << "Run07D";
                    } else if(modeltype[i] == "GRM"){
                        if(link[i] == "logit"){
                            d3Pidt2lmdu = zeros<mat>(npari, 2);
                            d4Pidt3klmdu = zeros<mat>(npari, 2);
                            d5hidt4duuniq = zeros<vec>(npari);
                            d4Pidt4klmn = zeros<mat>(p, 2);
                            arma::vec d3Pidt3klm = zeros<vec>(2);
                            arma::vec d2Pidt2lm = d2Pidt2.tube(templ, tempm);
                            //Define constants here
                            d3Pidt3klm(0) = aparsi(tempk) * (d2Pidt2(templ, tempm, 0) * (1.0 - 2.0 * Pi(0)) - 2.0 * dPidt(templ, 0) * dPidt(tempm, 0));
                            d3Pidt3klm(1) = aparsi(tempk) * (d2Pidt2(templ, tempm, 1) * (1.0 - 2.0 * Pi(1)) - 2.0 * dPidt(templ, 1) * dPidt(tempm, 1));
                           // Rcout << "Run07A";
                            d4hidt4uniq = aparsi(tempj) * (d3Pidt3klm(0) + d3Pidt3klm(1));
                            for(arma::uword u = 0; u < npari; u++){
                                d3Pidt2lmdu(u, 0) = aparsi(templ) * (d2Pidtdu(tempm, u, 0) * (1.0 - 2.0 * Pi(0)) - 2.0 * dPidt(tempm, 0) * dPidu(u, 0));
                                d3Pidt2lmdu(u, 1) = aparsi(templ) * (d2Pidtdu(tempm, u, 1) * (1.0 - 2.0 * Pi(1)) - 2.0 * dPidt(tempm, 1) * dPidu(u, 1));
                            }
                            d3Pidt2lmdu(templ, 0) += dPidt(tempm, 0) * (1.0 - 2.0 * Pi(0));
                            d3Pidt2lmdu(templ, 1) += dPidt(tempm, 1) * (1.0 - 2.0 * Pi(1));
                           // Rcout << "Run07B";
                            for(arma::uword u = 0; u < npari; u++){
                                d4Pidt3klmdu(u, 0) = aparsi(tempk) * (d3Pidt2lmdu(u, 0) * (1.0 - 2.0 * Pi(0)) - 2.0 * d2Pidt2(templ, tempm, 0) * dPidu(u, 0) - 2.0 * (d2Pidtdu(templ, u, 0) * dPidt(tempm, 0) + dPidt(templ, 0) * d2Pidtdu(tempm, u, 0)));
                                d4Pidt3klmdu(u, 1) = aparsi(tempk) * (d3Pidt2lmdu(u, 1) * (1.0 - 2.0 * Pi(1)) - 2.0 * d2Pidt2(templ, tempm, 1) * dPidu(u, 1) - 2.0 * (d2Pidtdu(templ, u, 1) * dPidt(tempm, 1) + dPidt(templ, 1) * d2Pidtdu(tempm, u, 1)));
                            }
                            d4Pidt3klmdu(tempk, 0) += d2Pidt2(templ, tempm, 0) * (1.0 - 2.0 * Pi(0)) - 2.0 * dPidt(templ, 0) * dPidt(tempm, 0);
                            d4Pidt3klmdu(tempk, 1) += d2Pidt2(templ, tempm, 1) * (1.0 - 2.0 * Pi(1)) - 2.0 * dPidt(templ, 1) * dPidt(tempm, 1);
                           // Rcout << "Run07C";
                            for(arma::uword u = 0; u < npari; u++){
                                d5hidt4duuniq(u) = aparsi(tempj) * (d4Pidt3klmdu(u, 0) + d4Pidt3klmdu(u, 1));
                            }
                            d5hidt4duuniq(tempj) += (d3Pidt3klm(0) + d3Pidt3klm(1));
                            for(arma::uword kk = 0; kk < p; kk++){
                                arma::uword tempn = 99;
                                for(arma::uword tt = 0; tt < pi; tt++){
                                    if(dimi(tt) == kk) tempn = tt; 
                                }
                                if(tempn == 99) continue;
                                arma::vec d3Pidt3lmn = zeros<vec>(2);
                                //Cross-check with d3gjdt3 (checks out)
                                d3Pidt3lmn(0) = aparsi(templ) * (d2Pidt2(tempm, tempn, 0) * (1.0 - 2.0 * Pi(0)) - 2.0 * dPidt(tempm, 0) * dPidt(tempn, 0));
                                d3Pidt3lmn(1) = aparsi(templ) * (d2Pidt2(tempm, tempn, 1) * (1.0 - 2.0 * Pi(1)) - 2.0 * dPidt(tempm, 1) * dPidt(tempn, 1));
                                
                                //Cross-check with d4gjdt4 (checks out)
                                d4Pidt4klmn(kk, 0) = aparsi(tempk) * (d3Pidt3lmn(0) * (1.0 - 2.0 * Pi(0)) - 2.0 * d2Pidt2(templ, tempm, 0) * dPidt(tempn, 0) - 2.0 * d2Pidt2(templ, tempn, 0) * dPidt(tempm, 0) - 2.0 * dPidt(templ, 0) * d2Pidt2(tempm, tempn, 0));
                                d4Pidt4klmn(kk, 1) = aparsi(tempk) * (d3Pidt3lmn(1) * (1.0 - 2.0 * Pi(1)) - 2.0 * d2Pidt2(templ, tempm, 1) * dPidt(tempn, 1) - 2.0 * d2Pidt2(templ, tempn, 1) * dPidt(tempm, 1) - 2.0 * dPidt(templ, 1) * d2Pidt2(tempm, tempn, 1));
                                d5hidt5uniq(kk) += aparsi(tempj) * (d4Pidt4klmn(kk, 0) + d4Pidt4klmn(kk, 1));
                            }
                        }
                        else if(link[i] == "probit"){
                            
                            d4hidt4uniq = d4hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) ;
                            for(arma::uword kk = 0; kk < p; kk++){
                                arma::uword tempn = 99;
                                for(arma::uword tt = 0; tt < pi; tt++){
                                    if(dimi(tt) == kk) tempn = tt; 
                                }
                                if(tempn == 99) continue;
                                d5hidt5uniq(kk) = d5hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * aparsi(tempn) ;
                            }
                            
                           // if( n==0 ){
                             //   Rcout << " uniq ready.";
                           // }
                            d5hidt4duuniq = zeros<vec>(npari);
                            // SJ: apar
                            for(arma::uword u = 0; u < pi; u++){
                                d5hidt4duuniq(u) = d5hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * thetai(u) + (aparsi(tempk) * aparsi(templ) * aparsi(tempm) + aparsi(tempj) * aparsi(templ) * aparsi(tempm) + aparsi(tempj) * aparsi(tempk) * aparsi(tempm) + aparsi(tempj) * aparsi(tempk) * aparsi(templ)) * d4hconst; // SJ: We need thetai to be the reduced dimension theta.
                            }
                            //if( n==0 ){
                             //   Rcout << " a ready.";
                            //}
                            // SJ: bpar
                            if(y(n, j) == 1){
                                d5hidt4duuniq(pi) += d5hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ;
                            }
                            else if(y(n, j) == mi(i)){
                                d5hidt4duuniq(npari - 1) += d5hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ; 
                            }
                            else{
                                d5hidt4duuniq(pi + y(n, j) - 2) += d5hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ;
                                d5hidt4duuniq(pi + y(n, j) - 1) += d5hconst * aparsi(tempj) * aparsi(tempk) * aparsi(templ) ;
                            }

                            //if( n==0 ){
                             //   Rcout << " b ready." << std::endl;
                            //}
                        }
                        //Rcout << "Run07D";
                    } else if(modeltype[i] == "negbin"){
						double topow = 1.0 + phi * explinpred;
                        d4hidt4uniq = (phi * ydouble + 1.0) * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm);
						//if(n == 0){
						//	Rcout << "d4hnd4i done" << std::endl;
						//}
						d5hidt4duuniq = zeros<vec>(npari);
						for(arma::uword u = 0; u < pi; u++){
                            double myeq1 = (tempj == u);
							double myeq2 = (tempk == u);
							double myeq3 = (templ == u);
							double myeq4 = (tempm == u);
                            d5hidt4duuniq(u) = (phi * ydouble + 1.0) * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * thetai(u) + (phi * ydouble + 1.0) * (explinpred * (2.0 * pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * thetai(u) - 4.0 * (phi * ydouble + 1.0) * (phi * exp2linpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 5.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * thetai(u) + myeq1 * (phi * ydouble + 1.0) * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 4.0) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) + myeq2 * (phi * ydouble + 1.0) * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(templ) * aparsi(tempm) + myeq3 * (phi * ydouble + 1.0) * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(tempm) + myeq4 * (phi * ydouble + 1.0) * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ);
                        }
						//if(n == 0){
						//	Rcout << "d5hnda done" << std::endl;
						//}
						d5hidt4duuniq(pi) = (phi * ydouble + 1.0) * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) + (phi * ydouble + 1.0) * (explinpred * (2.0 * pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) - 4.0 * (phi * ydouble + 1.0) * (phi * exp2linpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 5.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm);
						//if(n == 0){
						//	Rcout << "d5hndb done" << std::endl;
						//}
						d5hidt4duuniq(pi + 1) = ydouble * (explinpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) + (phi * ydouble + 1.0) * (explinpred * (2.0 * phi * exp2linpred - 4.0 * explinpred)) / pow(topow, 4.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) - 4.0 * (phi * ydouble + 1.0) * (exp2linpred * (pow(phi, 2.0) * exp2linpred - 4.0 * phi * explinpred + 1.0)) / pow(topow, 5.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm);
						//if(n == 0){
						//	Rcout << "d5hndphi done" << std::endl;
						//}
						for(arma::uword kk = 0; kk < p; kk++){
                            arma::uword tempn = 99;
                            for(arma::uword tt = 0; tt < pi; tt++){
                                if(dimi(tt) == kk) tempn = tt; 
                            }
                            if(tempn == 99) continue;
                            d5hidt5uniq(kk) += -(phi * ydouble + 1.0) * (explinpred * (pow(phi, 3.0) * exp3linpred - 11.0 * pow(phi, 2.0) * exp2linpred + 11.0 * phi * explinpred - 1.0)) / pow(topow, 5.0) * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * aparsi(tempn);
                        }
						//if(n == 0){
						//	Rcout << "d5hndu done" << std::endl;
						//}
					} else if(modeltype[i] == "poisson"){
						d4hidt4uniq = explinpred * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm);
						d5hidt4duuniq = zeros<vec>(npari);
						//if( n==0 ){
						//	Rcout << " d4hdz4 ready." << std::endl;
						//}
						for(arma::uword u = 0; u < pi; u++){
                            double myeq1 = (tempj == u);
							double myeq2 = (tempk == u);
							double myeq3 = (templ == u);
							double myeq4 = (tempm == u);
                            d5hidt4duuniq(u) = aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * explinpred * thetai(u) + explinpred * (aparsi(tempk) * aparsi(templ) * aparsi(tempm) * myeq1  + aparsi(tempj) * aparsi(templ) * aparsi(tempm) * myeq2 + aparsi(tempj) * aparsi(tempk) * aparsi(tempm) * myeq3 + aparsi(tempj) * aparsi(tempk) * aparsi(templ) * myeq4);
                        }
						//if( n==0 ){
						//	Rcout << " da ready." << std::endl;
						//}
						d5hidt4duuniq(pi) = aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * explinpred;
						//if( n==0 ){
						//	Rcout << " db ready." << std::endl;
						//}
                        for(arma::uword kk = 0; kk < p; kk++){
                            arma::uword tempn = 99;
                            for(arma::uword tt = 0; tt < pi; tt++){
                                if(dimi(tt) == kk) tempn = tt; 
                            }
                            if(tempn == 99) continue;
                            d5hidt5uniq(kk) += explinpred * aparsi(tempj) * aparsi(tempk) * aparsi(templ) * aparsi(tempm) * aparsi(tempn);
                        }
						//if( n==0 ){
						//	Rcout << " dz ready." << std::endl;
						//}
					} else if(modeltype[i] == "normal"){
						d5hidt4duuniq = zeros<vec>(npari);
					}
					//BA 2022-07-04: Add NRM
                    d4hdt4uniq(h, 0) += d4hidt4uniq;
                    for(arma::uword ii = 0; ii < p; ii++){
                        d4hdt4uniq(h, 1 + ii) += d5hidt5uniq(ii);
                    }
                     //if( n==0 ){
                     //    Rcout << " 1st assign ready.";
                     //}
                    //Derivatives with respect to item parameters
                    //Go from item i in group g to correct entry
                    //apars
					//BA 2022-06-29: Just copy derivatives for estimated pars?
					//BA 2022-07-02:
					for(arma::uword pari = 0; pari < npari; pari++){
						if(fulltounique(parindex + pari) == 0) continue;
						myjindex = fulltounique(parindex + pari) - 1;
						d4hdt4uniq(h, 1 + p + myjindex) += d5hidt4duuniq(pari);
					}
					
                     //if( n==0 ){
                     //    Rcout << " b assign ready.";
                     //}
                }
                 //if( n==0 ){
                 //    Rcout << " 4th order ready."  << std::endl;
                 //}
            }
            //Rcout << "Run08";
            //Need to define this per group, such that the starting value is from the first item parameter in group 'g'
            //Update is the same.
			//BA 2022-07-03: Updated 
            parindex += (G - group(n)) * npari;            
        }
        // SJ: End of item j.
		//BA: poisson and normal updated
        
        //if(n == 0){
        //Rcout << hn;
        //Rcout << dhn;
        //Rcout << d2hn;
        //Rcout << dhndu;
        //}
        //Two loops below is expensive part with many parameters, since matroperations? Can we simplify? How to time it..
        // SJ: I do not think there is a better way. We need the trace( matrix * matrix )
        arma::mat Bmat = inv(d2hn);
        arma::cube dBmatdu = zeros<cube>(p, p, nuniquepar);
        arma::vec tracedldu = zeros<vec>(nuniquepar);
        //Third derivatives w/ resp. to regression parameters are zero
        for(arma::uword u = 0; u < nuniquepar; u++){
            arma::mat tempmat11 = d3hndu(span(0, p - 1), span(0, p - 1), span(u, u));
            arma::mat tempmat22 = Bmat * tempmat11;
            arma::vec diagvec = tempmat22.diag();
            tracedldu(u) = sum(diagvec);
            dBmatdu(span(0, p - 1), span(0, p - 1), span(u, u)) = -tempmat22 * Bmat;
        }
        //Rcout << "Run09";
        arma::mat dtdu = -Bmat * d2hndu;
        arma::cube dBmatdt = zeros<cube>(p, p, p);
        arma::vec dldt = zeros<vec>(p);
        arma::vec tracedldt = zeros<vec>(p);
        for(arma::uword j = 0; j < p; j++){
            arma::mat tempmatA = d3hn(span(0, p - 1), span(0, p - 1), span(j));
            arma::mat tempmatB = Bmat * tempmatA;
            arma::vec diagvecB = tempmatB.diag();
            tracedldt(j) = sum(diagvecB);
            dBmatdt(span(0, p - 1), span(0, p - 1), span(j)) = -tempmatB * Bmat;
        }
        
        //First-order Laplace (hack solution)
        //double detd2hn;
        if(accuracy == 1){
            //detd2hn = det(d2hn);
			double logdetd2hn;
			double logdetd2hn_sign;
			log_det( logdetd2hn, logdetd2hn_sign, d2hn );
			//loglik(n) = -hn;
            loglik(n) = -hn - 0.5 * logdetd2hn + 0.5 * pdouble * log2pi;
            dldt = -dhn - 0.5 * tracedldt;
            gloglik.col(n) = -dhndu - 0.5 * tracedldu + dtdu.t() * dldt;
            // if ( n == 0 ){
            //     //Rcout << "dtdu = " << dtdu << std::endl;
            //     Rcout << "tracedldt = " << tracedldt.t() << std::endl;
            //     Rcout << "dhn = " << dhn.t() << std::endl;
            //     Rcout << "dldt = " << dldt.t() << std::endl;
            // }
        }
        //Second-order Laplace
		//BA 2022-06-29: Below is a general implementation, not affected by fixed pars or different item types
        if( accuracy == 2){
            //Rcout << "Run10";
            double R21 = 0.0;
            double R22 = 0.0;
            double R23 = 0.0;
            arma::vec depsilonfgdu = zeros<vec>(nuniquepar);
            arma::vec depsilonfgdt = zeros<vec>(p); 
            for(arma::uword h = 0; h < nUniqComb_3rd_4; h++){
                arma::uword iposjkm = UniqComb_3rd_4(h, 0);
                arma::uword iposrst = UniqComb_3rd_4(h, 1);
                arma::uword j = UniqComb_3rd_4(h, 2);
                arma::uword k = UniqComb_3rd_4(h, 3);
                arma::uword l = UniqComb_3rd_4(h, 4);
                arma::uword r = UniqComb_3rd_4(h, 5);
                arma::uword s = UniqComb_3rd_4(h, 6);
                arma::uword ti = UniqComb_3rd_4(h, 7);
                
                R23 += UniqComb_3rd_4(h, 8) * d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 0) * Bmat(j, k) * Bmat(l, r) * Bmat(s, ti);
                //Gradient
                for(arma::uword u = 0; u < nuniquepar - nregpar; u++) depsilonfgdu(u) += UniqComb_3rd_4(h, 8) * (1.0 / 8.0) * ((d3hdt3uniq(iposjkm, 1 + p + u) * d3hdt3uniq(iposrst, 0) + d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 1 + p + u)) * Bmat(j, k) * Bmat(l, r) * Bmat(s, ti) + d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 0) * (dBmatdu(j, k, u) * Bmat(l, r) * Bmat(s, ti) + Bmat(j, k) * dBmatdu(l, r, u) * Bmat(s, ti) + Bmat(j, k) * Bmat(l, r) * dBmatdu(s, ti, u)));
                for(arma::uword ii = 0; ii < p; ii++) depsilonfgdt(ii) += UniqComb_3rd_4(h, 8) * (1.0 / 8.0) * ((d3hdt3uniq(iposjkm, 1 + ii) * d3hdt3uniq(iposrst, 0) + d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 1 + ii)) * Bmat(j, k) * Bmat(l, r) * Bmat(s, ti) + d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 0) * (dBmatdt(j, k, ii) * Bmat(l, r) * Bmat(s, ti) + Bmat(j, k) * dBmatdt(l, r, ii) * Bmat(s, ti) + Bmat(j, k) * Bmat(l, r) * dBmatdt(s, ti, ii)));
            }
            //Rcout << "Run11";
            for(arma::uword h = 0; h < nUniqComb_3rd_6; h++){
                arma::uword iposjkm = UniqComb_3rd_6(h, 0);
                arma::uword iposrst = UniqComb_3rd_6(h, 1);
                arma::uword j = UniqComb_3rd_6(h, 2);
                arma::uword k = UniqComb_3rd_6(h, 3);
                arma::uword l = UniqComb_3rd_6(h, 4);
                arma::uword r = UniqComb_3rd_6(h, 5);
                arma::uword s = UniqComb_3rd_6(h, 6);
                arma::uword ti = UniqComb_3rd_6(h, 7);
                R22 += UniqComb_3rd_6(h, 8) * d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 0) * Bmat(j, k) * Bmat(l, r) * Bmat(s, ti);
                for(arma::uword u = 0; u < nuniquepar - nregpar; u++) depsilonfgdu(u) += UniqComb_3rd_6(h, 8) * (1.0 / 12.0) * ((d3hdt3uniq(iposjkm, 1 + p + u) * d3hdt3uniq(iposrst, 0) + d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 1 + p + u)) * Bmat(j, k) * Bmat(l, r) * Bmat(s, ti) + d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 0) * (dBmatdu(j, k, u) * Bmat(l, r) * Bmat(s, ti) + Bmat(j, k) * dBmatdu(l, r, u) * Bmat(s, ti) + Bmat(j, k) * Bmat(l, r) * dBmatdu(s, ti, u)));
                for(arma::uword ii = 0; ii < p; ii++) depsilonfgdt(ii) += UniqComb_3rd_6(h, 8) * (1.0 / 12.0) * ((d3hdt3uniq(iposjkm, 1 + ii) * d3hdt3uniq(iposrst, 0) + d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 1 + ii)) * Bmat(j, k) * Bmat(l, r) * Bmat(s, ti) + d3hdt3uniq(iposjkm, 0) * d3hdt3uniq(iposrst, 0) * (dBmatdt(j, k, ii) * Bmat(l, r) * Bmat(s, ti) + Bmat(j, k) * dBmatdt(l, r, ii) * Bmat(s, ti) + Bmat(j, k) * Bmat(l, r) * dBmatdt(s, ti, ii)));
            }
            //Rcout << "Run12";
            for(arma::uword h = 0; h < nUniqComb; h++){
                arma::uword ipos = UniqComb(h, 0);
                arma::uword j = UniqComb(h, 1);
                arma::uword k = UniqComb(h, 2);
                arma::uword l = UniqComb(h, 3);
                arma::uword m = UniqComb(h, 4);
                
                R21 += UniqComb(h, 5) * d4hdt4uniq(ipos, 0) * Bmat(j, k) * Bmat(l, m);
                for(arma::uword ii = 0; ii < p; ii++) depsilonfgdt(ii) += -UniqComb(h, 5) * 0.125 * (d4hdt4uniq(ipos, 1 + ii) * Bmat(j, k) * Bmat(l, m) + d4hdt4uniq(ipos, 0) * (dBmatdt(j, k, ii) * Bmat(l, m) + Bmat(j, k) * dBmatdt(l, m, ii)));
                for(arma::uword u = 0; u < nuniquepar - nregpar; u++) depsilonfgdu(u) += -UniqComb(h, 5) * 0.125 * (d4hdt4uniq(ipos, 1 + p + u) * Bmat(j, k) * Bmat(l, m) + d4hdt4uniq(ipos, 0) * (dBmatdu(j, k, u) * Bmat(l, m) + Bmat(j, k) * dBmatdu(l, m, u)));
            }	
            //Rcout << "Run13";
            epsilon = 0.0;
            epsilon = -(1.0 / 2.0) * ((1.0 / 4.0) * R21 - ((1.0 / 6.0) * R22 + (1.0 / 4.0) * R23));
            gepsilon = depsilonfgdu / (1.0 + epsilon);
            tepsilon = depsilonfgdt / (1.0 + epsilon);
			double logdetd2hn;
			double logdetd2hn_sign;
			log_det( logdetd2hn, logdetd2hn_sign, d2hn );
            //detd2hn = det(d2hn);
            loglik(n) = -hn - 0.5 * logdetd2hn + 0.5 * pdouble * log2pi + log(1.0 + epsilon);
            dldt = -dhn - 0.5 * tracedldt + tepsilon;
            gloglik.col(n) = -dhndu - 0.5 * tracedldu + gepsilon + dtdu.t() * dldt;
        }
    }
    //Rcout << "RunXX";
    return Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                              Rcpp::Named("gradient") = gloglik);
    // //Rcpp::Named("map") = mapmatr,
    // //Rcpp::Named("nc") = failind);
}

//[[Rcpp::export]]
arma::mat chol_deri ( arma::mat mat_in ) {
    //  SJ: Function to compute the derivative of the Cholesky decomposition.
    int row_no = mat_in.n_rows ;
    mat output = mat_in ;
    for( int i=0; i<row_no; i++ ){
        for( int j=i; j<row_no; j++ ){
            if( i==j ){
                output(i,j) = 0.5*mat_in(i,j) ;
            }
            else {
                output(i,j) = 0.0 ;
            }
        }
    }
    return(output);
}



// [[Rcpp::export]]
List mglogLGrad_quad(arma::vec pars, arma::mat estfixpars, arma::mat y, arma::mat theta, arma::uword J, arma::vec mi, arma::uword p, Rcpp::List model, std::vector<std::string> modeltype, std::vector<std::string> link, arma::uword N, arma::mat covstruct, arma::uword G, arma::vec group, Rcpp::List filters, arma::mat X, bool adapt, bool fullexp, arma::mat quadp, arma::vec quadw, std::string method, arma::uvec npartype) {
    
    // SJ: The function inputs are mostly the same as the inputs of mglogLGrad, expect the new ones
    //     bool adapt: whether it is adaptive quadrature or not
    //     bool fullexp: whether it is fully exponential or not. true corresponds to direct maximization.
    //                  We actually want to approximate integral / integral.
    //                  The possible combinations of adapt and fullexp are:
    //                    adapt = false && fullexp = false: regular approximation to both integrals
    //                    adapt = true && fullexp = false: adaptive approximation to both integrals
    //                    adapt = false && fullexp = true: not supported! Results do not mean anything.
    //                    adapt = true && fullexp = true: fully exponential approximation to the whole ratio.
    //     arma::mat quadp: a matrix of quadrature points. Here quadrature can be GH quadrature or (quasi) monte carlo.
    //     arma::vec quadw: a vector of quadrature weights.
    //     std::string method: which approximation method is used, ghq or mc.
    //     The plan is that
    //       regular Gauss-Hermite quadrature if adapt = false && method = "ghq"
    //       adaptive Gauss-Hermite quadrature if adapt = true && method = "ghq"
    //       regular monte carlo if adapt = false && method = "mc"
    //       adaptive monte carlo if adapt = true && method = "mc"
    //     To determine whether it is regular monte carlo or quasi-monte carlo, we should prepare quadrature weights properly in R.
    //       C++ does not distinguish different monte carlo methods. It only uses the supplied quadrature points and weights. 
    if( adapt == false && fullexp == true ){
        Rcout << "adapt = false && fullexp = true is currently not supported. Switch to adapt = true && fullexp = true instead." << std::endl ;
        adapt = true ;
    }
    arma::uword noQuad = static_cast<uword>(quadw.size()); // SJ: Number of quadrature points.
	
	//BA 2022-07-08: Updated, consistent with mglogLGrad
	arma::uword nitpar = npartype(0);
	arma::uword nregpar = npartype(1);
	arma::uword nmupar = npartype(2);
	arma::uword nvarpar = npartype(3);
	arma::uword ncovpar = npartype(4);
	arma::uword nuniquepar = nitpar + nregpar + nmupar + nvarpar + ncovpar;
	
    Rcpp::List apars(G);
    Rcpp::List bpars(G);
    arma::uword parindex = 0;
	arma::uvec fulltounique = conv_to<uvec>::from(estfixpars.col(2));
	
	Rcpp::List mylists = tabletolist(estfixpars, J, mi, model, modeltype, G);
	apars = mylists(0);
	bpars = mylists(1);
	arma::vec tmpobj = mylists(2);
	arma::uvec npar = conv_to<uvec>::from(tmpobj);
	
	arma::uword nitpartot = G * sum(npar);
	arma::uword nregpartot = X.n_cols - 1;
	nregpartot = G * nregpartot * p;
	arma::uword nmupartot = G * p;
	arma::uword nvarpartot = G * p;
	arma::uword ncovpartot = G * p * (p - 1) / 2;
	
	//BA 2022-07-03: This is the total number of "active" parameters (the number of rows in estfixpars)
	arma::uword npartot = nitpartot + nregpartot + nmupartot + nvarpartot + ncovpartot;
	
	//Initialize means and covariance matrices
    arma::mat mu = zeros<mat>(p, G);
    arma::cube Sigma = zeros<cube>(p, p, G);
	
	//Rcout << "Run001";
	//BA 2022-07-02: We just run through all the entries
	//BA 2022-07-05: We update the index here to the point where the distribution parameters begin
	parindex = G * sum(npar);	
    for(arma::uword g = 0; g < G; g++){
		mu.col(g) = estfixpars(span(parindex, parindex + p - 1), span(1));
        parindex += p;			
    }
	//Rcout << "Run002";
	for(arma::uword g = 0; g < G; g++){
		Sigma.slice(g) = diagmat(estfixpars(span(parindex, parindex + p - 1), span(1)));
        parindex += p;
	}
	for(arma::uword g = 0; g < G; g++){
        //fill in covariance parameters
		for(arma::uword j = 0; j < (p - 1); j++){
			for(arma::uword k = j + 1; k < p; k++){
				Sigma(j, k, g) = estfixpars(parindex, 1);
				Sigma(k, j, g) = estfixpars(parindex, 1);
				parindex += 1;
			}
        }
    }

	/*
    //Define filters.
    //Input is the object 'filters', a list with entries:
    //uniqi3, uniqi4, Uniq3, Uniq4, UniqComb_3rd_4, UniqComb_3rd_6, UniqComb (0, 1, 2, 3, 4, 5, 6)
    //List w/ J entries, List w/ J entries, matrix, matrix, matrix, matrix, matrix 
    Rcpp::List uniqi3 = filters(0); // SJ: unique 3rd-order derivatives that we need to compute.
    Rcpp::List uniqi4 = filters(1); // SJ: unique 4th-order derivatives that we need to compute.
    arma::mat Uniq3 = filters(2);
    arma::mat Uniq4 = filters(3);
    arma::mat UniqComb_3rd_4 = filters(4);
    arma::mat UniqComb_3rd_6 = filters(5);
    arma::mat UniqComb = filters(6);
    
    arma::uword nUniqComb_3rd_4 = static_cast<uword>(UniqComb_3rd_4.n_rows);
    arma::uword nUniqComb_3rd_6 = static_cast<uword>(UniqComb_3rd_6.n_rows);
    arma::uword nUniqComb = static_cast<uword>(UniqComb.n_rows);
	*/
	//Number of regression parameters
    arma::uword nbetapar = X.n_cols - 1;
    nbetapar = nbetapar * p;
			
    arma::mat betamat;	
    arma::mat mydhdb;
    arma::cube myd2hdtdb;
	//BA 2022-07-08: This is currently broken (wrong indices)
    if(G == 1){
        if(nregpartot > 0){
            betamat = zeros<mat>(p, nregpartot / p + 1);
            arma::uword betaparind = npartot - nregpartot;
            for(arma::uword tt = 0; tt < p; tt++){
                betamat(span(tt, tt), span(1, nregpartot / p)) = trans(pars(span(betaparind, betaparind + nregpartot / p - 1)));
                betaparind += nregpartot / p;
            }
            //Rcout << "Run001";
            for(arma::uword ss = 0; ss < p; ss++){
                arma::mat temp1 = X(span(0, N - 1), span(1, nregpartot / p));
                arma::vec temp2 = trans(betamat(span(ss, ss), span(1, nregpartot / p)));
                betamat(span(ss, ss), 0) = -mean(temp1 * temp2);
            }
            //Rcout << "Run002";
            mydhdb = dhdb(theta, Sigma.slice(0), betamat, X, p, N);
            //Rcout << "Run003";
            myd2hdtdb = d2hdtdb(theta, Sigma.slice(0), betamat, X, p, N);
            //Rcout << "Run004";
        }
    }
    
    //Rcout << "Run02";
    //Define inverses of covariance matrices
    arma::cube invSigma = zeros<cube>(p,p,G);
    for(arma::uword g = 0; g < G; g++) invSigma.slice(g) = inv(Sigma.slice(g));
    //Define constants
    arma::vec logsqrtdet2pisigma = zeros<vec>(G);
    for(arma::uword g = 0; g < G; g++){
        arma::mat tmpsig = 2.0 * M_PI * Sigma.slice(g);
        double templogsqrtdet2pisigma = det(tmpsig);
        templogsqrtdet2pisigma = sqrt(templogsqrtdet2pisigma);
        templogsqrtdet2pisigma = log(templogsqrtdet2pisigma);
        logsqrtdet2pisigma(g) = templogsqrtdet2pisigma;
    }
	
    //These are the output variables.
    //Output log-likelihood and gradient for each individual.
    arma::vec loglik = zeros<vec>(N);
    arma::mat gloglik = zeros<mat>(nuniquepar, N);
    
    arma::vec addvec1;
    arma::vec addvec2;
	arma::uword muindex;
	arma::uword varindex;
	arma::uword covindex;
    for(arma::uword n = 0; n < N; n++){
		//BA 2022-07-05: Enable escaping from C++ to R prompt.
		Rcpp::checkUserInterrupt();
        //Prepare all the things we need.
        //For each item, we add to the entries of these objects.
        arma::uword ng = group(n);
        arma::vec mung = mu.col(ng);
        arma::mat invSigmang = invSigma.slice(ng);
        double logsqrtdet2pisigmag = logsqrtdet2pisigma(ng);
        //double hn = 0.0;
        //arma::vec dhn = zeros<vec>(p);
        arma::mat d2hn = zeros<mat>(p, p);
        arma::cube d3hn = zeros<cube>(p, p, p);
        //arma::vec dhndu = zeros<vec>(nuniquepar);
        arma::mat d2hndu = zeros<mat>(p, nuniquepar);
        arma::cube d3hndu = zeros<cube>(p, p, nuniquepar);
        if(nregpartot > 0) mung = betamat * trans(X.row(n));
        arma::vec temptheta = trans(theta.row(n));
        arma::vec tempvec = temptheta - mung;
		//Rcout << "Run003";
		if(nmupar > 0){
            //arma::vec dhdmu = zeros<vec>(p);
            arma::mat d2hdtdmu = zeros<mat>(p, p);
			for(arma::uword g = 0; g < G; g++){
				if(group(n) != g) continue;
				//dhdmu(span(0, p - 1)) = -invSigma.slice(g) * tempvec;
				d2hdtdmu(span(0, p - 1), span(0, p - 1)) = -invSigma.slice(g);
			}
			//BA 2022-07-03: We loop over all mu-parameters in estfixpars, but only add estimated ones to the entries in dhndu and d2hndu
			for(arma::uword u = 0; u < p; u++){
				if(fulltounique(nitpartot + group(n) * p + u) == 0) continue;
				muindex = fulltounique(nitpartot + group(n) * p + u) - 1;
				//dhndu(span(muindex, muindex)) += dhdmu(u);
				d2hndu(span(0, p - 1), span(muindex, muindex)) += d2hdtdmu.col(u);
			}
        }
		
		//Rcout << "Run004";
        //arma::vec dhds = zeros<vec>(p + p * (p - 1) / 2 );
        arma::mat d2hdtds = zeros<mat>(p, p + p * (p - 1) / 2);
        arma::cube d3hdt2ds = zeros<cube>(p, p, p + p * (p - 1) / 2);
        
        //Definition of checkmat1 and checkmat2 from covarianceprep()
        //arma::mat checkmat1 = covstruct;
		//BA 2022-07-04: We always go through all entries - will waste some time but makes code easier. covstruct input not needed, but can keep for now to avoid issues elsewhere.
        arma::mat checkmat2 = zeros<mat>(p, 2);
        checkmat2 = join_cols(checkmat2, covstruct);
        for(arma::uword j = 0; j < p; j++){
            checkmat2(j, 0) = j;
            checkmat2(j, 1) = j;
        }
        arma::mat checkmat;
        arma::uword s;
        arma::mat ident = eye(p, p);
        //Derivatives for covariance matrix parameters
		//Rcout << "Run005";
		//BA 2022-07-02: This needs some housekeeping - but should work essentially as-is
		//BA 2022-07-03: Fixed, is it OK?
		for(arma::uword g = 0; g < G; g++){
            if(group(n) != g) continue;
            checkmat = checkmat2;
            arma::mat tempmat = -0.5 * (2.0 * invSigma.slice(g) - invSigma.slice(g) % ident -  2.0 * invSigma.slice(g) * tempvec * tempvec.t() * invSigma.slice(g) + invSigma.slice(g) * tempvec * tempvec.t() * invSigma.slice(g) % ident);
			arma::mat tempmat1MG;
            //Rcout << "Run101";			
			//Variance parameters
			s = 0;
			for(arma::uword u = 0; u < p; u++){
				arma::uword u1 = checkmat(u, 0);
                arma::uword u2 = checkmat(u, 1);
                tempmat1MG = zeros<mat>(p, p);
                //dhds(s) = -tempmat(u1, u2);
                tempmat1MG(u1, u2) = 1.0;
                d2hdtds.col(s) = -invSigma.slice(g) * tempmat1MG * invSigma.slice(g) * tempvec;
                s += 1;
			}
			//Covariance parameters
			for(arma::uword u = p; u < p + p * (p - 1) / 2; u++){
				arma::uword u1 = checkmat(u, 0);
                arma::uword u2 = checkmat(u, 1);
				tempmat1MG = zeros<mat>(p, p);
				//dhds(s) = -tempmat(u1, u2);
				tempmat1MG(u1, u2) = 1.0;
                tempmat1MG(u2, u1) = 1.0;
                d2hdtds.col(s) = -invSigma.slice(g) * tempmat1MG * invSigma.slice(g) * tempvec;
                s += 1;
            }
        }
	
		//New code
		//Rcout << "Run006";
		//BA 2022-07-02: This needs some housekeeping - but should work essentially as-is
		for(arma::uword g = 0; g < G; g++){
            if(group(n) != g) continue;
            checkmat = checkmat2;
            arma::vec tempvec = temptheta - mung;
			s = 0;
			for(arma::uword u = 0; u < p; u++){
				arma::uword u1 = checkmat(u, 0);
                arma::uword u2 = checkmat(u, 1);
                arma::mat tempmatMG2 = zeros<mat>(p, p);
				tempmatMG2(u1, u2) = 1.0;
                d3hdt2ds(span(0, p - 1), span(0, p - 1), span(s, s)) = - invSigma.slice(g) * tempmatMG2 * invSigma.slice(g);
                s += 1;
			}
			for(arma::uword u = p; u < p + p * (p - 1) / 2; u++){
				arma::uword u1 = checkmat(u, 0);
                arma::uword u2 = checkmat(u, 1);
                arma::mat tempmatMG2 = zeros<mat>(p, p);
				tempmatMG2(u1, u2) = 1.0;
                tempmatMG2(u2, u1) = 1.0;
                d3hdt2ds(span(0, p - 1), span(0, p - 1), span(s, s)) = - invSigma.slice(g) * tempmatMG2 * invSigma.slice(g);
                s += 1;
			}
		}
		
       // Rcout << "Run04";
		//BA 2022-07-04: Separate placement of variance/covariance derivatives since they are ordered separately by group
		//BA 2022-07-03: Link variance parameters between estfixpars and unique pars
		//BA 2022-07-03: Identify correct parameter from estfixpars, and link to unique parameter (if estimated)
		for(arma::uword u = 0; u < p; u++){
			if(fulltounique(nitpartot + nmupartot + group(n) * p + u) == 0) continue;
			varindex = fulltounique(nitpartot + nmupartot + group(n) * p + u) - 1;
			//dhndu(span(varindex, varindex)) += dhds(u);
            d2hndu(span(0, p - 1), span(varindex, varindex)) += d2hdtds.col(u);
            d3hndu(span(0, p - 1), span(0, p - 1), span(varindex, varindex)) += d3hdt2ds.slice(u);
		}
		//if( n == 0 ){
		//	Rcout << "Index: " << varindex << std::endl;
        //}
		//BA 2022-07-03: Link covariance parameters between estfixpars and unique pars
		//BA 2022-07-03: Identify correct parameter from estfixpars, and link to unique parameter (if estimated)
		for(arma::uword u = p; u < p + p * (p - 1) / 2; u++){
			if(fulltounique(nitpartot + nmupartot + nvarpartot + group(n) * p * (p - 1) / 2 + u - p) == 0) continue;
			covindex = fulltounique(nitpartot + nmupartot + nvarpartot + group(n) * p * (p - 1) / 2 + u - p) - 1;
			//dhndu(span(covindex, covindex)) += dhds(u);
            d2hndu(span(0, p - 1), span(covindex, covindex)) += d2hdtds.col(u);
            d3hndu(span(0, p - 1), span(0, p - 1), span(covindex, covindex)) += d3hdt2ds.slice(u);
		}
		//if( n == 0 ){
		//	Rcout << "Index: " << covindex << std::endl;
        //}
	
        //Derivatives for regression parameters (only for single group case right now)
		//BA 2022-07-04: This needs updating, broken order now
        if(nregpartot != 0){
            //dhndu(span(nitpar, nitpar + nregpar - 1)) = trans(mydhdb.row(n));
            arma::mat myd2hdtdbn = myd2hdtdb(span(n, n), span(0, p - 1), span(0, nregpar - 1));
            d2hndu(span(0, p - 1), span(nitpar, nregpar - 1)) = -myd2hdtdbn;
        }
		
        //Define object sizes for unique third- and fourth-order derivatives
        //For each item, we will add the unique entries to these objects
        //These will then be combined in the final calculation of the likelihood
        //The order is:  third- or fourth-order derivative, derivatives with respect to the latent variables, derivatives with respect to the unknown parameters
		/*
        arma::uword nuniq3all = static_cast<uword>(Uniq3.n_rows);
        arma::uword nuniq4all = static_cast<uword>(Uniq4.n_rows);
        arma::uword ncolsuniq34 = 1 + p + nuniquepar;
        arma::mat d3hdt3uniq = zeros<mat>(nuniq3all, ncolsuniq34);
        arma::mat d4hdt4uniq = zeros<mat>(nuniq4all, ncolsuniq34);
        */
		
        //Multivariate normal distribution
        arma::mat tmpp = -0.5 * trans(temptheta - mung) * invSigmang * (temptheta - mung);
        //hn += -(tmpp(0,0) - logsqrtdet2pisigmag);
        //dhn += trans(invSigmang) * (temptheta - mung);
        d2hn += invSigmang;
		//BA 2022-07-02: Here, we need to keep track of the new order of things in the object estfixpars
		//BA 2022-07-03: We start with first item, which are the first parameters in estfixpars
        parindex = 0;
        arma::uword myjindex;
        Rcpp::List aparsg = apars(ng);
        Rcpp::List bparsg = bpars(ng);
        //--SJ: If we are using GRM with probit link, we need to compute a constant term for third order derivative. 
        double d3hconst = 0.0;
        
		//BA 2022-07-08: Updated up to here (see line 2071 in lamle.cpp)
        //Add unique stuff to objects that are later combined and weighted
        //Everything in this loop depends on the different IRT models
        // if( n==0 ){
        //     Rcout << "Start measurement model" ;
        // }
		//BA 2023-05-16: Start adding negbin, normal, and poisson
		//BA 2023-05-16: For AGHQ, we add to d2hn, d3hn, d2hndu, d3hndu
        for(arma::uword i = 0; i < J; i++){
            double midouble = mi(i);
			double sumPimi = 0.0;
			double ydouble = y(n, i);
			double phi = 1.0;
            arma::vec mival = regspace(1.0, midouble);
            //  SJ: p = dimension of latent variables, including those with zero slopes.
            //      pj = number of parameters in a: reduced dimension.
            arma::uvec dimi = model(i);
            arma::uword pi = static_cast<uword>(dimi.n_elem);
            arma::uword npari = npar(i);
			parindex += group(n) * npari;
            
            //Missing value handling, just skip to the next item. Missing coded as 9999.
            if(y(n, i) == 9999){
                parindex += (G - group(n)) * npari;
                continue;
            }
            //dimi defines the non-zero 3rd order derivatives for item i
            arma::vec aparsi = aparsg(i);
            arma::vec bparsi = bparsg(i);
            arma::vec thetai = temptheta(dimi);
			double linpred = bparsi(1) + sum(aparsi % thetai);
			double toexp2 = 2.0 * linpred;
			double explinpred = exp(linpred);
			double exp2linpred = exp(toexp2);
            arma::vec Pi;
            //  SJ: OBS! For GRM with probit link, we are computing logPi, not Pi.
			//	BA: Also, output is different between GRM and GPCM? Should be, for efficiency. Derivatives for GPCM and NRM have same structure, however.
			// BA: Objects below not needed for some models
            arma::mat dPidt;
            arma::cube d2Pidt2;
            arma::mat dPidu;
            arma::cube d2Pidtdu;
            // SJ: Added for new function (item-wise)
			// BA: Only for GRM with probit link
            arma::mat dlogPrdt;
            arma::mat d2logPrdt2;
            arma::cube d3logPrdt3;
            arma::mat dlogPrdu;
            arma::mat d2logPrdtdu;
            arma::cube d3logPrdt2du;
            //Add NRM to the functions below
            //  SJ: Compute the gradient of h with respect to item parameters

            if(modeltype[i] == "GPCM"){
                Pi = gi(thetai, aparsi, bparsi, modeltype[i], mi(i), ydouble);
                dPidt = dgidz(thetai, aparsi, bparsi, modeltype[i], Pi, mi(i), pi, ydouble);
                d2Pidt2 = d2gidz2(thetai, aparsi, bparsi, modeltype[i], Pi, dPidt, mi(i), pi, ydouble);
                dPidu = dgidu(thetai, aparsi, bparsi, modeltype[i], Pi, mi(i), pi, ydouble);
                d2Pidtdu = d2gidzdu(thetai, aparsi, bparsi, modeltype[i], Pi, dPidt, dPidu, mi(i), pi, ydouble);
                sumPimi = sum(Pi % mival);
            }
            else if(modeltype[i] == "GRM"){
                // SJ: I have not yet optimize the code. It contains lots of repeated and unnecessary calculations!
				// BA: OK! It is indeed much slower than my previous code (for GRM).
                if(link[i] == "logit"){
                    //  SJ: Using arma::mat is expected to be faster. I am simply using Rcpp::List for simplicity.
                    //Rcpp::List GRM_List = dgjdu_GRM_string(y(n, jj), thetajj, aparsjj, bparsjj, link[jj], mj(jj), pj, nparj) ;
                    Rcpp::List GRM_List = item_GRM(y(n, i), thetai, aparsi, bparsi, link[i], mi(i), pi, npari) ;
                    //arma::vec Pi_temp = GRM_List["Pi"];
                    //Pi = Pi_temp;
                    //Pi = gj_GRM_logit(y(n, jj), thetajj, aparsjj, bparsjj, mj(jj)) ;
                    arma::mat dPidt_temp = GRM_List["dPidt"];
                    dPidt = dPidt_temp; // SJ: I need dPidt.
                    arma::cube d2Pidt2_temp = GRM_List["d2Pidt2"];
                    d2Pidt2 = d2Pidt2_temp;
                    arma::mat dPidu_temp = GRM_List["dPidu"];
                    dPidu = dPidu_temp;
                    arma::cube d2Pidtdu_temp = GRM_List["d2Pidtdu"];
                    d2Pidtdu = d2Pidtdu_temp;
                    
                    // SJ: I have used one function for 1st and 2nd order derivatives, and one function for 3rd order derivatives.
                    arma::mat GRMmat = item_GRM_quad(y(n, i), thetai, aparsi, bparsi, link[i], mi(i), pi, npari) ;
                    Pi = zeros<vec>(2);
                    Pi(0) = GRMmat(0,0);
                    Pi(1) = GRMmat(1,0);
                    //hn += -log(Pi(0) - Pi(1));
                    dlogPrdt = GRMmat.submat(0, 1, pi - 1, 1);
                    dlogPrdu = GRMmat.submat(pi, 1, pi + npari -1, 1);
                    //dhnduj = -1.0 * dlogPrdu;
                    arma::cube GRM3rdcube = item_GRM_3rd(y(n, i), thetai, aparsi, bparsi, link[i], mi(i), pi, npari, Pi, true) ;
                    d3logPrdt3 = GRM3rdcube( span(0, pi - 1), span(0, pi - 1), span(0, pi - 1) );
                    d3logPrdt2du = GRM3rdcube( span(0, pi - 1), span(0, pi - 1), span(pi, pi + npari - 1) );
                    d2logPrdt2 = GRMmat.submat(0, 2, pi -1, 1 + pi);
                    d2logPrdtdu = GRMmat.submat(pi, 2, pi + npari -1, 1 + pi).t();
                    
                }
                else if(link[i] == "probit"){
                    Pi = gj_GRM_probit(y(n, i), thetai, aparsi, bparsi, mi(i)) ;//  SJ: Not only Pr().
                    //   = Pivec( span(0,1) ) ;
                    dPidu = dgjdu_GRM_probit(y(n, i), thetai, Pi(4), Pi(5), mi(i), pi, npari);
                    dPidt = dgjdt_GRM_probit(aparsi, Pi(4), Pi(5), pi);
                    d2Pidt2 = d2gjd2t_GRM_probit(y(n, i), aparsi, Pi(2), Pi(3), Pi(4), Pi(5), mi(i), pi);
                    d2Pidtdu = d2gjdtdu_GRM_probit(y(n, i), thetai, aparsi, Pi(2), Pi(3), Pi(4), Pi(5), mi(i), pi, npari);
                    d3hconst = d3gjd3t_GRM_probit( y(n, i), aparsi, Pi(2), Pi(3), Pi(4), Pi(5), mi(i), pi );
                    //hn += -log(Pi(0) - Pi(1));
                    //dhnduj = -1.0 * (dPidu.col(0) - dPidu.col(1));
                }
            }
            
            //  //  SJ: Copy the gradient for item j into the gradient of all parameters
            //  //      This part does not depend on the item type.
            //  //apars
            //  for(arma::uword parj = 0; parj < pj; parj++){
            //      myjindex = fulltounique(parindex + parj);
            //      dhndu(myjindex) += dhnduj(parj);
            //  }
            //  //bpars
            //  for(arma::uword parj = pj; parj < (pj + mj(jj) - 1); parj++){
            //      myjindex = fulltounique(parindex + pj * (mj(jj) - 1) + parj - pj);
            //      dhndu(myjindex) += dhnduj(parj);
            //  }
            
            //  SJ: The outer for-loop compute dh / dtheta (dhn)
            //  SJ: Why not compute the following loop in external modeltype specific functions?
            //      For the item parameters, we need dlogPr / dtheta
            //Rcout << "Run05";
            //Index variables
            //  SJ: If I understand the code correctly, index1 to index3 are used to find out which column of the complete apar we are talking about.
            //      Sometimes, we are only working with the nonzero entries in apar(jj).
            //      But here, we need to work with the entire apar(jj) vector.
            arma::uword index1 = 0;
            arma::uword index2;
            arma::uword index3;
            //This loop does repeated calculations for multidimensional models and can be improved. However, we need to fill out these objects in the end so changing it might not be worth it.
            for(auto j : dimi){
                index2 = 0;
                //posvec is defined as 1.0 for the value such that index1 is equal to the parameter value, i.e. the index1-th entry 
                arma::vec posvec = zeros<vec>(npari);
                addvec1 = zeros<vec>(npari);
                //-- SJ: if(modeltype[jj] == 2 || modeltype[jj] == 4){
                if(modeltype[i] == "GPCM"){
                    //dhn(i) += -aparsjj(index1) * (ydouble - sumPimj);
                    posvec(index1) = 1.0;
                    for(arma::uword v = 0; v < npari; v++) addvec1(v) = aparsi(index1) * sum(trans(dPidu.row(v)) % mival);
                    addvec1 += -posvec * (ydouble - sumPimi);
                } else if(modeltype[i] == "GRM"){
                    if(link[i] == "logit"){
                        addvec1 = conv_to< arma::vec >::from( -1.0 * d2logPrdtdu.row(index1) ); 
                    }
                    else if(link[i] == "probit"){
                        for(arma::uword v = 0; v < npari; v++) addvec1(v) = aparsi(index1) * (dPidu(v, 0) + dPidu(v, 1));
                        addvec1(index1) += -(Pi(4) - Pi(5));
                    }
                } else if(modeltype[i] == "negbin"){
					phi = bparsi(2);
					for(arma::uword v = 0; v < pi; v++){
						double myeq = (v == index1);
						addvec1(v) = -myeq * ydouble + (phi * ydouble + 1.0) * (explinpred) / (1.0 + phi * explinpred) * aparsi(index1) * thetai(v) - (phi * ydouble + 1.0) * (phi * exp2linpred) / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(index1) * thetai(v) + myeq * (phi * ydouble + 1.0) * explinpred / (1.0 + phi * explinpred) * (aparsi(index1) / aparsi(v));
					}
					addvec1(pi) = (phi * ydouble + 1.0) * explinpred / (1.0 + phi * explinpred) * aparsi(index1) - (phi * ydouble + 1.0) * (phi * exp2linpred) / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(index1);
					addvec1(pi + 1) = ydouble * (explinpred) / (1.0 + phi * explinpred) * aparsi(index1) - (phi * ydouble + 1.0) * exp2linpred / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(index1);
					//if(n == 0){
					//	Rcout << "d2hndu done" << std::endl;
					//}
				} else if(modeltype[i] == "normal"){
					phi = bparsi(2);
					for(arma::uword v = 0; v < pi; v++){
						double myeq = (v == index1);
						addvec1(v) = (aparsi(index1) * thetai(v) + (linpred - ydouble) * myeq) / phi;
					}
					addvec1(pi) = aparsi(index1) / phi;
					//addvec1(pi + 1) = -aparsi(index1) * (linpred - ydouble) / (phi * phi);
					addvec1(pi + 1) = ydouble * aparsi(index1) / pow(phi, 2.0) - linpred * aparsi(index1) / pow(phi, 2.0);
				} else if(modeltype[i] == "poisson"){
					for(arma::uword v = 0; v < pi; v++){
						double myeq = (index1 == v);
						addvec1(v) = aparsi(index1) * explinpred * thetai(v)  + (explinpred - ydouble) * myeq;
					}
					addvec1(pi) = aparsi(index1) * explinpred;
				}
                
                for(arma::uword pari = 0; pari < npari; pari++){
					if(fulltounique(parindex + pari) == 0) continue;
					myjindex = fulltounique(parindex + pari) - 1;
					d2hndu(j, myjindex) += addvec1(pari);
				}
                //Rcout << "Run24";
                
                //  SJ: The inner for-loop computes d2h / d2theta (d2hn), the Hessian of h with respect to latent variables
                for(auto k : dimi){
                    index3 = 0;
                    addvec1 = zeros<vec>(npari);
                    //-- SJ: if(modeltype[jj] == 2 || modeltype[jj] == 4){
                     if(modeltype[i] == "GPCM"){
                        d2hn(j, k) += aparsi(index1) * sum(trans(dPidt.row(index2)) % mival);
                        arma::mat mytemp1 = d2Pidtdu(span(index2), span(0, npari - 1), span(0, mi(i) - 1));
                        for(arma::uword v = 0; v < npari; v++) addvec1(v) = aparsi(index1) * sum(trans(mytemp1.row(v)) % mival);
                        addvec1 += posvec * sum(trans(dPidt.row(index2)) % mival);
                    } else if(modeltype[i] == "GRM"){
                        if(link[i] == "logit"){
                            d2hn(j, k) += -d2logPrdt2(index1, index2);
                            arma::rowvec direct_subset = d3logPrdt2du(arma::span(index1,index1), arma::span(index2,index2), arma::span::all);
                            addvec1 = conv_to< arma::vec >::from( -1.0 * direct_subset );
                        }
                        else if(link[i] == "probit"){
                            d2hn(j, k) += -d2Pidt2(index1,index2,0);
                            //  SJ:  How to make this addvec1 right??
							//  BA:  Hm, this is probit - not sure about what's going on. Is it OK?
                            //for(arma::uword v = 0; v < nparj; v++) addvec1(v) = d2Pidtdu(index2, v, 0);
                            //addvec1(index1) += -d2Pidt2(i,j,0);
                        }
                    } else if(modeltype[i] == "negbin"){
						double topow  = 1.0 + phi * exp(linpred);
						d2hn(j, k) += (phi * ydouble + 1.0) * (explinpred) / pow(topow, 2.0) * aparsi(index1) * aparsi(index2);
						for(arma::uword v = 0; v < pi; v++){
							double myeq1 = (index1 == v);
							double myeq2 = (index2 == v);
							addvec1(v) = (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * aparsi(index1) * aparsi(index2) * thetai(v) - 2.0 * (phi * ydouble + 1.0) * (phi * exp2linpred) / pow(topow, 3.0) * aparsi(index1) * aparsi(index2) * thetai(v) + myeq1 * (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * aparsi(index2) + myeq2 * (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * aparsi(index1);
						}
						addvec1(pi) = (phi * ydouble + 1.0) * explinpred / pow(topow, 2.0) * aparsi(index1) * aparsi(index2) - 2.0 * (phi * ydouble + 1.0) * (phi * exp2linpred) / pow(topow, 3.0) * aparsi(index1) * aparsi(index2);
						addvec1(pi + 1) = ydouble * explinpred / pow(topow, 2.0) * aparsi(index1) * aparsi(index2) - 2.0 * (phi * ydouble + 1.0) * exp2linpred / pow(topow, 3.0) * aparsi(index1) * aparsi(index2);
						//if(n == 0){
						//	Rcout << "d3hndu done" << std::endl;
						//}
					} else if(modeltype[i] == "normal"){
						d2hn(j, k) += aparsi(index1) * aparsi(index2) / phi;
						for(arma::uword v = 0; v < pi; v++){
							double myeq1 = (index1 == v);
							double myeq2 = (index2 == v);
							addvec1(v) = (aparsi(index2) * myeq1 + aparsi(index1) * myeq2) / phi;
						}
						//Derivative with respect to the intercept is zero
						addvec1(pi + 1) = -aparsi(index1) * aparsi(index2) / pow(phi, 2.0);
					} else if(modeltype[i] == "poisson"){
						d2hn(j, k) += explinpred * aparsi(index1) * aparsi(index2);
						for(arma::uword v = 0; v < pi; v++){
							double myeq1 = (index1 == v);
							double myeq2 = (index2 == v);
							addvec1(v) = aparsi(index1) * aparsi(index2) * explinpred * thetai(v) + explinpred * (aparsi(index2) * myeq1 + aparsi(index1) * myeq2);
						}
						addvec1(pi) = aparsi(index1) * aparsi(index2) * explinpred;
					}
                    for(arma::uword pari = 0; pari < npari; pari++){
						if(fulltounique(parindex + pari) == 0) continue;
						myjindex = fulltounique(parindex + pari) - 1;
						d3hndu(j, k, myjindex) += addvec1(pari);
					}
                    
                    // SJ: the last for-loop to compute d3h / d3theta (d3hn)
                    for(auto l : dimi){
                        //-- SJ: if(modeltype[jj] == 2 || modeltype[jj] == 4){
                        if(modeltype[i] == "GPCM"){
                            arma::vec mytemp2 = d2Pidt2(span(index2), span(index3), span(0, mi(i) - 1));
                            d3hn(j, k, l) += aparsi(index1) * sum(mytemp2 % mival);
                        } else if(modeltype[i] == "GRM"){
                            if(link[i] == "logit"){
                                //d3hn(i, j, k) += aparsjj(index1) * (d2Pidt2(index2, index3, 0) + d2Pidt2(index2, index3, 1));
                                d3hn(j, k, l) += -d3logPrdt3(index1, index2, index3);
                            }
                            else if(link[i] == "probit"){
                                d3hn(j, k, l) += d3hconst * aparsi(index1) * aparsi(index2) * aparsi(index3) ;
                            }
                        } else if(modeltype[i] == "negbin"){
							d3hn(j, k, l) += -(phi * ydouble + 1.0) * explinpred * (phi * explinpred - 1.0) / ((1.0 + phi * explinpred) * (1.0 + phi * explinpred) * (1.0 + phi * explinpred)) * aparsi(index1) * aparsi(index2) * aparsi(index3); 
						} else if(modeltype[i] == "poisson"){
							d3hn(j, k, l) += explinpred * aparsi(index1) * aparsi(index2) * aparsi(index3);
						}
                        index3 += 1;
                    }
                    index2 += 1;
                }
                index1 += 1;
            }
            
            //Need to define this per group, such that the starting value is from the first item parameter in group 'g'
            //BA 2022-07-03: Updated 
            parindex += (G - group(n)) * npari;
        }
        //if( n==0 ){
        //    Rcout << " Finish measurement model";
        //}

        //Two loops below is expensive part with many parameters, since matroperations? Can we simplify? How to time it..
        // SJ: I do not think there is a better way. We need the trace( matrix * matrix )
        arma::mat Bmat = inv(d2hn);
        // SJ: If adaptive quadrature, we need to work with the cholesky decomposition.
        arma::mat lmat = zeros<mat>(p, p) ; 
        arma::cube dcholdu;
        arma::cube dcholdt;
        if( adapt == true ){
            // SJ: Cholesky decomposition, organized as a lower triangular matrix.
            lmat = arma::chol(Bmat, "lower") ;
            // SJ: Both dcholdu and dcholdt are only needed in the fully exponential approximation.
            if( fullexp == true ){
                //  Compute the third-order derivative, which is needed to accounted for the dependence between the factor scores and the parameters. 
                dcholdu = zeros<cube>(p, p, nuniquepar);
                dcholdt = zeros<cube>( p, p, p ) ;
                //  Here, compute dsqrt(2) * L / dtheta. This is needed to account for the dependence between the adaptive quadrature points and the mode, in turn the parameters.
                for( arma::uword j = 0; j < p; j++){
                    dcholdt.slice(j) = sqrt(2.0) * lmat * chol_deri( -1.0 * lmat.t() * d3hn.slice(j) * lmat ) ;
                }
            }
        }
        //if( n==0 ){
        //    Rcout << " dcholdt OK" << std::endl;
        //}
        // SJ: Compute the derivative of the negative Hessian matrix.
        arma::cube dBmatdu = zeros<cube>(p, p, nuniquepar);
        arma::vec tracedldu = zeros<vec>(nuniquepar);
        //Third derivatives w/ resp. to regression parameters are zero, if we observe all regressors in the latent regression
        for(arma::uword u = 0; u < nuniquepar; u++){
            arma::mat tempmat11 = d3hndu(span(0, p - 1), span(0, p - 1), span(u, u));
            arma::mat tempmat22 = Bmat * tempmat11;
            arma::vec diagvec = tempmat22.diag();
            tracedldu(u) = sum(diagvec); // SJ: why not use trace()? BA: Probably because I didn't think of it..
            dBmatdu(span(0, p - 1), span(0, p - 1), span(u, u)) = -tempmat22 * Bmat;
            
            // SJ: If we are working with adaptive quadrature, we need to get the derivative of the Cholesky decomposition with respect to the parameters.
            if( adapt == true && fullexp == true ){
                dcholdu.slice(u) = sqrt(2.0) * trans( lmat * chol_deri( -1.0 * lmat.t() * tempmat11 * lmat ) ) ;
            }
        }
        //if( n==0 ){
        //    Rcout << " dcholdu OK" << std::endl;
        //}
        //Rcout << "Run09";
        arma::mat dtdu = -Bmat * d2hndu;
        //arma::cube dBmatdt = zeros<cube>(p, p, p); // SJ: Do we really need this guy? BA: No, we do not need for AGHQ.
        arma::vec dldt = zeros<vec>(p);
        arma::vec tracedldt = zeros<vec>(p);
        for(arma::uword j = 0; j < p; j++){
            arma::mat tempmatA = d3hn(span(0, p - 1), span(0, p - 1), span(j)); // SJ: why not d3hn.slice(i)? BA: Some issue that I can't remember, but maybe works fine now.
            arma::mat tempmatB = Bmat * tempmatA;
            arma::vec diagvecB = tempmatB.diag();
            tracedldt(j) = sum(diagvecB);
            //dBmatdt(span(0, p - 1), span(0, p - 1), span(i)) = -tempmatB * Bmat;
        }
        
        // ----------------------------------------------------------------------------- //
        // ----------------------------------------------------------------------------- //
        // SJ: Start to compute the likelihood function and the gradient.
        //Rcout << "Start quadrature loop" << std::endl ;
        double sumwexph = 0.0 ; // SJ: w * exp(h)
        arma::mat wdexphdt = zeros<mat>(p, 1);
        arma::mat dadaptdtheta; // SJ: This is needed if adapt == true && fullexp == true
        arma::vec dhn_quad ; // SJ: This is dh / dtheta, needed only if adapt == true && fullexp == true 
        arma::vec quad_grad = zeros<vec>(nuniquepar);
        for( arma::uword qi = 0; qi < noQuad; qi++ ){
            
            // SJ: The current quadrature points, depending on bool adapt.
            arma::mat quadp_nquad = quadp.col(qi); 
            // SJ: If adapt = false, we only need to extract each column in quadp.
            arma::mat theta_nquad = quadp_nquad;
            // SJ: Otherwise, we need to translate and dilate quadrature points.
            if( adapt == true ){
                theta_nquad = sqrt(2.0) * lmat * quadp_nquad + temptheta; // SJ: temptheta should be the mode.
            }
            // SJ: Initiate derivative
            arma::vec dhndu_quad = zeros<vec>(nuniquepar); 
            if( adapt == true && fullexp == true ){
                dhn_quad = zeros<vec>(p);
            }
            
            // SJ: Multivariate normal distribution for h. 
            //     OBS: Only valid for single group.
			// BA: Why is it only valid for single group?
			// BA: Covariance pars OK? Mupars are missing? Just add mupar code and is OK? (implemented)
            arma::mat tempvec_quad = theta_nquad - mung;
            arma::mat tmpp_quad = -0.5 * tempvec_quad.t() * invSigmang * tempvec_quad; // SJ: use as_scalar?
            double hn_quad = -(tmpp_quad(0,0) - logsqrtdet2pisigmag);
            //BA: New addition
			//Need change here for some means/vars estimated 
			/*
			if(nmupar > 0){
				arma::vec dhdmu = zeros<vec>(nmupar); // BA: This is dh / d(mean parameter)
				if(G != 1){
					for(arma::uword g = 1; g < G; g++){
						if(group(n) != g) continue;
						arma::uword s = (g - 1) * p;
						dhdmu(span(s, s + p - 1)) = -invSigma.slice(g) * tempvec_quad;
					}
				}
				dhndu_quad(span(nitpar, nitpar + nmupar - 1)) = dhdmu;			
			}
			*/
			
			if(nmupar > 0){
				arma::vec dhdmu = zeros<vec>(p);
				for(arma::uword g = 0; g < G; g++){
					if(group(n) != g) continue;
					dhdmu(span(0, p - 1)) = -invSigma.slice(g) * tempvec_quad;
				}
				for(arma::uword u = 0; u < p; u++){
					if(fulltounique(nitpartot + group(n) * p + u) == 0) continue;
					muindex = fulltounique(nitpartot + group(n) * p + u) - 1;
					dhndu_quad(span(muindex, muindex)) += dhdmu(u);
				}
			}
            // SJ: Derivatives for covariance matrix parameters
            // SJ: Do not fully understand the code for covariance parameters.
            //arma::vec dhds_quad = zeros<vec>(nsigmapar); // SJ: This is dh / d(covariance parameter)
			arma::vec dhds_quad = zeros<vec>(p + p * (p - 1) / 2 );
			arma::mat checkmat2 = zeros<mat>(p, 2);
			checkmat2 = join_cols(checkmat2, covstruct);
			for(arma::uword j = 0; j < p; j++){
				checkmat2(j, 0) = j;
				checkmat2(j, 1) = j;
			}
			arma::mat checkmat;
			arma::uword s;
			arma::mat ident = eye(p, p);
			//Need change here for some means/vars estimated 
            for(arma::uword g = 0; g < G; g++){
                if(group(n) != g) continue;
                // SJ: Do not understand checkmat1 and checkmat2. Do we need to re-assign checkmat1 and checkmat2? BA: I don't remember exactly why I did it like this, but it works.. I think there was some trial-and-error involved here which made it convoluted. Probably a bit wasteful to do it like this but I suspect it does not contribute much to runtime
                //if(g == 0) checkmat = checkmat1; else checkmat = checkmat2;
				checkmat = checkmat2;
                arma::mat tempmat = -0.5 * (2.0 * invSigma.slice(g) - invSigma.slice(g) % ident -  2.0 * invSigma.slice(g) * tempvec_quad * tempvec_quad.t() * invSigma.slice(g) + invSigma.slice(g) * tempvec_quad * tempvec_quad.t() * invSigma.slice(g) % ident);
                //Rcout << "Run101";
                //s denotes which index to start from
                //start from number of variance parameters of previous groups plus the number of covariance parameters of previous groups
                //for g = 1, start from number of rows in covstruct (equal to number of covariance parameters)
                //for g = 2, start from number of rows in covstruct times two plus p
                //for g = 3, start from number of rows in covstruct times three plus p times two
                //Second derivatives variance parameters?
                //Add check for no covariance parameters?
                //Equality constraints across groups?
                //if(g == 0) s = 0; else s = (g - 1) * p + g * ncovgpar;
                //Variance parameters
				s = 0;
				for(arma::uword u = 0; u < p; u++){
					arma::uword u1 = checkmat(u, 0);
					arma::uword u2 = checkmat(u, 1);
					dhds_quad(s) = -tempmat(u1, u2);
					s += 1;
				}
				//Covariance parameters
				for(arma::uword u = p; u < p + p * (p - 1) / 2; u++){
					arma::uword u1 = checkmat(u, 0);
					arma::uword u2 = checkmat(u, 1);
					dhds_quad(s) = -tempmat(u1, u2);
					s += 1;
				}
            }
            // SJ: If we have free parameters in the covariance matrix of latent variables, ...
            for(arma::uword u = 0; u < p; u++){
				if(fulltounique(nitpartot + nmupartot + group(n) * p + u) == 0) continue;
				varindex = fulltounique(nitpartot + nmupartot + group(n) * p + u) - 1;
				dhndu_quad(span(varindex, varindex)) += dhds_quad(u);
			}
			//if( n == 0 ){
			//	Rcout << "Index: " << varindex << std::endl;
			//}
			//BA 2022-07-03: Link covariance parameters between estfixpars and unique pars
			//BA 2022-07-03: Identify correct parameter from estfixpars, and link to unique parameter (if estimated)
			for(arma::uword u = p; u < p + p * (p - 1) / 2; u++){
				if(fulltounique(nitpartot + nmupartot + nvarpartot + group(n) * p * (p - 1) / 2 + u - p) == 0) continue;
				covindex = fulltounique(nitpartot + nmupartot + nvarpartot + group(n) * p * (p - 1) / 2 + u - p) - 1;
				dhndu_quad(span(covindex, covindex)) += dhds_quad(u);
			}
			//if( n == 0 ){
			//	Rcout << "Index: " << covindex << std::endl;
			//}
            // SJ: dh / d(latent variable) due to the normal distribution assumption.
            if( adapt == true && fullexp == true ){
                dhn_quad = invSigmang * tempvec_quad;
                //if ( n == 0 ){
                //    Rcout << "cov dhn_quad = " << dhn_quad << std::endl;
                //}
            }
            
            //Derivatives for regression parameters (only for single group case right now)
            // SJ: This is not programmed yet!!!!!!!!!
			// BA: Yeah, we need to add support for regression w/ MG (wasn't implemented in my previous code either, for MG, need to define restrictions in a different way with possible implications for gradient computation)
			//BA 2022-07-08: Need to fix regression parameters, new index
            if(nregpartot != 0){
                dhndu_quad(span(nitpar + nregpartot, nitpar + nregpartot - 1)) = trans(mydhdb.row(n));
            }
            
            // SJ: Measurement model
            // SJ: I think I need a new parindex here, which is used the as the starting point.
            arma::uword parindex_quad = 0;
            for(arma::uword i = 0; i < J; i++){
				double midouble = mi(i);
				double ydouble = y(n, i);
				double phi = 1.0;
				arma::vec mival = regspace(1.0, midouble);
				//  SJ: p = dimension of latent variables, including those with zero slopes.
				//      pi = number of latent variables with non-zero slopes
				arma::uvec dimi = model(i);
				arma::uword pi = static_cast<uword>(dimi.n_elem);
				//BA 2022-07-02: Need to be item-specific
				//BA 2022-07-03: Updated, we have a new object "npari" which gives the number of *non-zero* parameters for each item (same number in each group)
				arma::uword npari = npar(i);
				parindex_quad += group(n) * npari;
			    //if( n == 0 && qi == 0){
				//    Rcout << "item " << i << " pi=" << pi << " mi(i)=" << mi(i) << " npari=" << npari << std::endl;
				//	Rcout << "Index: " << parindex_quad << std::endl;
				//}
				//Missing value handling, just skip to the next item. Missing coded as 0.
				//BA 2022-07-02: Need to change below based on group, since we have ordered the item parameters by item by group
				//BA 2022-07-03: Should be fixed now
				if(y(n, i) == 9999){
					parindex_quad += (G - group(n)) * npari;
					continue;
				}
				//dimi defines the non-zero 3rd order derivatives for item i
				arma::vec aparsi = aparsg(i);
				arma::vec bparsi = bparsg(i);
				arma::vec thetai = theta_nquad(dimi);
				arma::vec Pi;
				//Is the problem the linpred objects here? Should I use a different theta?
				double linpred = bparsi(1) + sum(aparsi % thetai);
				double explinpred = exp(linpred);
				/*
                double mjdouble = mj(jj);
                arma::vec mjval = regspace(1.0, mjdouble);
                //  SJ: p = dimension of latent variables, including those with zero slopes.
                //      pj = number of parameters in a
                //Change line below for non-GPCM or non-GRM items.
                arma::uvec dimj = model(jj);
                arma::uword pj = static_cast<uword>(dimj.n_elem);
                arma::uword nparj = pj + mj(jj) - 1;
                arma::uword nparjj = pj * (mj(jj) - 1) + (mj(jj) - 1);
                //Missing value handling, just skip to the next item. Missing coded as 0.
                if(y(n, jj) == 0){
					//BA: Corrected.
                    //parindex += nparjj;
					parindex_quad += nparjj;
                    continue;
                }
                //dimj defines the non-zero 3rd order derivatives for item jj
                arma::vec aparsjj = aparsg(jj);
                arma::vec bparsjj = bparsg(jj);
                arma::vec thetajj = theta_nquad(dimj); // SJ: I think we should redefine thetajj.
                arma::vec Pi;
				
				*/
                //  SJ: OBS! For GRM with probit link, we are computing logPi, not Pi.
                arma::mat dPidt;
                arma::mat dPidu;
                // SJ: Added for new function (item-wise)
                arma::mat dlogPrdt;
                arma::mat dlogPrdu;
                //Add NRM to the functions below
                //  SJ: Compute the gradient of h with respect to item parameters
                arma::vec dhndui_quad;
                if(modeltype[i] == "GPCM"){
                    Pi = gi(thetai, aparsi, bparsi, modeltype[i], mi(i), ydouble);
                    dPidt = dgidz(thetai, aparsi, bparsi, modeltype[i], Pi, mi(i), pi, ydouble);
                    dPidu = dgidu(thetai, aparsi, bparsi, modeltype[i], Pi, mi(i), pi, ydouble);
                    
                    double tempobj = Pi(y(n, i) - 1);
                    //ydouble = y(n, jj);
                    //sumPimj = sum(Pi % mjval);
                    hn_quad += -log(tempobj);
                    //Define range for entries into du-objects based on the placement of parameters for item jj
                    dhndui_quad = -dPidu.col(y(n, i) - 1) / tempobj;
                    
                    // SJ: If adapt == true && fullexp == true, we also need to compute dlog(Pr) / dtheta.
                    //     Is it the right way to specify the gradient with respect to theta?? <- Hm..?
                    if( adapt == true && fullexp == true ){
                        arma::mat dlogPrdt_quad = dPidt.col(y(n, i) - 1) / tempobj;
                        for(arma::uword ui = 0; ui < pi; ui++){
                            dhn_quad(dimi(ui), 0) += -dlogPrdt_quad(ui, 0);
                        }
                        //if( n == 0 ){
                        //    Rcout << "    item " << jj << " dhn part = " << -dlogPrdt_quad << std::endl;
                        //}
                    }
                    
                } else if(modeltype[i] == "GRM"){
                    if(link[i] == "logit"){
                        //  SJ: Using arma::mat is expected to be faster. I am simply using Rcpp::List for simplicity.
                        //Rcpp::List GRM_List = dgjdu_GRM_string(y(n, jj), thetajj, aparsjj, bparsjj, link[jj], mj(jj), pj, nparj) ;
                        //Rcpp::List GRM_List = itemj_GRM(y(n, jj), thetajj, aparsjj, bparsjj, link[jj], mj(jj), pj, nparj) ;
                        //arma::vec Pi_temp = GRM_List["Pi"];
                        //Pi = Pi_temp;
                        //Pi = gj_GRM_logit(y(n, jj), thetajj, aparsjj, bparsjj, mj(jj)) ;
                        //arma::mat dPidt_temp = GRM_List["dPidt"];
                        //dPidt = dPidt_temp; // SJ: I need dPidt.
                        //arma::mat dPidu_temp = GRM_List["dPidu"];
                        //dPidu = dPidu_temp;
                        
                        // SJ: I have used one function for 1st and 2nd order derivatives, and one function for 3rd order derivatives.
                        arma::mat GRMmat = item_GRM_quad(y(n, i), thetai, aparsi, bparsi, link[i], mi(i), pi, npari) ;
                        Pi = zeros<vec>(2);
                        Pi(0) = GRMmat(0,0);
                        Pi(1) = GRMmat(1,0);
                        hn_quad += -log(Pi(0) - Pi(1));
                        dlogPrdu = GRMmat.submat(pi, 1, pi + npari -1, 1);
                        dhndui_quad = -1.0 * dlogPrdu;
                        
                        // SJ: If adapt == true && fullexp == true, we also need to compute dlog(Pr) / dtheta.
                        if( adapt == true && fullexp == true ){
                            arma::mat dlogPrdt_quad = GRMmat.submat(0, 1, pi - 1, 1);
                            for(arma::uword ui = 0; ui < pi; ui++){
                                dhn_quad(dimi(ui), 0) += -dlogPrdt_quad(ui, 0);
                            }
                            //if( n == 0 ){
                            //    Rcout << "    item " << jj << " dhn part = " << -dlogPrdt_quad << std::endl;
                            //}
                        }
                    } else if(link[i] == "probit"){
                        Pi = gj_GRM_probit(y(n, i), thetai, aparsi, bparsi, mi(i)) ;//  SJ: Not only Pr(). BA: What's the output here?
                        //   = Pivec( span(0,1) ) ;
                        dPidu = dgjdu_GRM_probit(y(n, i), thetai, Pi(4), Pi(5), mi(i), pi, npari);
                        dPidt = dgjdt_GRM_probit(aparsi, Pi(4), Pi(5), pi); // SJ: dgjdt_GRM_probit actually computes dlogPi / dtheta, but to the lower dimension.
                        hn_quad += -log(Pi(0) - Pi(1));
                        dhndui_quad = -1.0 * (dPidu.col(0) - dPidu.col(1));
                        
                        // SJ: If adapt == true && fullexp == true, we also need to compute dlog(Pr) / dtheta.
                        if( adapt == true && fullexp == true ){
                            arma::mat dlogPrdt_quad = dPidt.col(0) - dPidt.col(1);
                            for(arma::uword ui = 0; ui < pi; ui++){
                                dhn_quad(dimi(ui), 0) += -dlogPrdt_quad(ui, 0);
                            }
                        }
                        
                    }
                } else if(modeltype[i] == "negbin"){
					phi = bparsi(2);
					double negbintolog1 = 1.0 / phi + explinpred;
					double negbintolog2 = 1.0 + phi * explinpred;
					double todigam1 = ydouble + 1.0 / phi;
					double todigam2 = 1.0 / phi;
					double togammaf = ydouble + 1.0;
					hn_quad += -lgamma(todigam1) + lgamma(togammaf) + lgamma(todigam2) - ydouble * linpred + ydouble * log(negbintolog1) + (1.0 / phi) * log(negbintolog2);
					dhndui_quad = zeros<vec>(npari);
					for(arma::uword pari = 0; pari < pi; pari++){
						dhndui_quad(pari) = -ydouble * thetai(pari) + (ydouble + 1.0 / phi) * (phi * explinpred) / (1.0 + phi * explinpred) * thetai(pari);
					}
					dhndui_quad(pi) = -ydouble + (ydouble + 1.0 / phi) * (phi * explinpred) / (1.0 + phi * explinpred);
					dhndui_quad(pi + 1) = -log(negbintolog2) / (phi * phi) + (ydouble + 1.0 / phi) * (explinpred) / (1.0 + phi * explinpred) - ydouble / phi + R::digamma(todigam1) / (phi * phi) - R::digamma(todigam2) / (phi * phi);
					
					if(adapt == true && fullexp == true){
						for(arma::uword ui = 0; ui < pi; ui++){
							dhn_quad(dimi(ui), 0) += -ydouble * aparsi(ui) + (phi * ydouble + 1.0) * explinpred / (1.0 + phi * explinpred) * aparsi(ui);
						}
					}
				} else if(modeltype[i] == "normal"){
					phi = bparsi(2);
					hn_quad += -(ydouble * linpred - linpred * linpred / 2.0 ) / phi + ydouble * ydouble / (2.0 * phi) + (log2pi + log(phi)) / 2.0;
					dhndui_quad = zeros<vec>(npari);
					for(arma::uword pari = 0; pari < pi; pari++){
						dhndui_quad(pari) = thetai(pari) * (linpred - ydouble) / phi;
					}
					dhndui_quad(pi) = (linpred - ydouble) / phi;
					dhndui_quad(pi + 1) = ydouble * linpred / pow(phi, 2.0) - pow(linpred, 2.0) / (2.0 * pow(phi, 2.0)) - pow(ydouble, 2.0) / (2.0 * pow(phi, 2.0)) + 1.0 / (2.0 * phi);
					
					if(adapt == true && fullexp == true){
						for(arma::uword ui = 0; ui < pi; ui++){
							dhn_quad(dimi(ui), 0) += aparsi(ui) * (linpred - ydouble) / phi;
						}
					}
				} else if(modeltype[i] == "poisson"){
					double yfac = rcpp_factorial(ydouble);
					hn_quad += -log(pow(explinpred, ydouble) * exp(-explinpred) / yfac);
					dhndui_quad = zeros<vec>(npari);
					for(arma::uword pari = 0; pari < pi; pari++){
						dhndui_quad(pari) = thetai(pari) * (explinpred - ydouble);
					}
					dhndui_quad(pi) = explinpred - ydouble;	
					if(adapt == true && fullexp == true){
						for(arma::uword ui = 0; ui < pi; ui++){
							dhn_quad(dimi(ui), 0) += (explinpred - ydouble) * aparsi(ui);
						}
					}
				}
        
                //  SJ: Copy the gradient for item j into the gradient of all parameters
                //      This part does not depend on the item type.
				for(arma::uword pari = 0; pari < npari; pari++){
				//BA 2022-07-02: Need to account for some parameters not estimated
					if(fulltounique(parindex_quad + pari) == 0) continue;
					myjindex = fulltounique(parindex_quad + pari) - 1;
					dhndu_quad(myjindex) += dhndui_quad(pari);
				}
				
                //Rcout << "Run08";
                //Need to define this per group, such that the starting value is from the first item parameter in group 'g'
                //Update is the same.
                parindex_quad += (G - group(n)) * npari;   
                
            }
            
            // SJ: Accumulate the quadrature formula.
            double wexph = quadw(qi) * exp(-1.0 * hn_quad) ;
            sumwexph += wexph ; 

            // SJ: In order to compute the gradient, the estimated mode depends on the parameters, 
            //     if adapt = true and fullexp = true.
            if( adapt == true && fullexp == true ){
                dadaptdtheta = eye<mat>(p, p) ;
                for ( arma::uword j = 0; j < p; j++) {
                    // SJ: We consider the derivative of the cholesky decomposition outside of the loop, including the coefficient sqrt(2).
                    dadaptdtheta.col(j) += dcholdt.slice(j) * quadp_nquad ;
                }
                wdexphdt += wexph * ( -1.0 * dadaptdtheta.t() * dhn_quad ) ;
                
                // SJ: We need to add the derivative because of the cholesky decomposition.
                //     In the current model, do we need to consider beta?
				// Yes, but my previous code didn't support combination of MG and regression-pars.
                //for(arma::uword u = 0; u < ntotpar - nbetapar; u++){
                for(arma::uword u = 0; u < nuniquepar; u++){
                    // SJ: We probably need a submat() after slice.
                    dhndu_quad(u) += as_scalar( quadp_nquad.t() * dcholdu.slice(u) * dhn_quad ) ;
                }
            }
            
            //if ( n == 0 ){
            //    Rcout << "dhn_quad = " << dhn_quad << std::endl;
            //}
            
            // SJ: The gradient approximated by quadrature.
            quad_grad += -1.0 * wexph * dhndu_quad ;
            
        }
        // SJ: Finish the likelihood function and the gradient.
        //Rcout << "Finish quadrature loop" << std::endl ;
        // ----------------------------------------------------------------------------- //
        // ----------------------------------------------------------------------------- //
        
        // SJ: log-likelihood value
        loglik(n) += log(sumwexph);
        if( adapt == true ){
            // SJ: log_det is more precise than log(det()). Maybe it does not matter since we do not work with large matrices. BA: Sounds good, I wasn't aware of this function, can change elsewhere too. 
            double logdetd2hn;
            double logdetd2hn_sign;
			//BA: Should be log_det of lmat ?
			//log_det( logdetd2hn, logdetd2hn_sign, lmat );
            log_det( logdetd2hn, logdetd2hn_sign, d2hn );
            loglik(n) += 0.5 * p * log(2.0) - 0.5 * logdetd2hn;
        }
        
        // SJ: gradient with respect to all parameters.
        //     If adapt == true && fullexp == true, we need to compute dlog(L) / dtheta 
        if( adapt == false && fullexp == false ){
            gloglik.col(n) = quad_grad / sumwexph;
        }
        if( adapt == true && fullexp == false ){
            gloglik.col(n) = quad_grad / sumwexph;
        }
        if( adapt == true && fullexp == true ){
            // SJ: we need to compute dldt, the gradient of approximated logL with respect to latent variables.
            dldt = wdexphdt / sumwexph - 0.5 * tracedldt ; // SJ: Does dimension/object type match? BA: Well, you suggest that it works? :)
            gloglik.col(n) = quad_grad / sumwexph - 0.5 * tracedldu + dtdu.t() * dldt ;
            //if ( n == 0 ){
            //    //Rcout << "dtdu = " << dtdu << std::endl;
            //    Rcout << "dhn = " << wdexphdt / sumwexph << std::endl;
            //    Rcout << "tracedldt = " << tracedldt.t() << std::endl;
            //    Rcout << "dldt = " << dldt.t() << std::endl;
            //}
        }
        
    }
    //Rcout << "RunXX";
    return Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                              Rcpp::Named("gradient") = gloglik);   
}
