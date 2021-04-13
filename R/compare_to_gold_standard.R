# This is a function for comparing an objects hash value to its gold standard.
# I am putting it here since some tests still don't use the modular testing function.
# Eventually I will move the function to the test directory.

compare_to_gold_standard = function(check_name,obj_to_hash,test_directory = 'C:/Test_DACOMP/'){
  hash_computation_result = digest::digest(obj_to_hash, algo="md5")
  gold_standard_dir = paste0(test_directory,'Gold_Standard_Values/')
  if(!dir.exists(gold_standard_dir)){
    dir.create(gold_standard_dir)
  }
  file_to_check = paste0(gold_standard_dir,check_name,'.RData')
  if(file.exists(file_to_check)){
    load(file = file_to_check) #=> hash_gold_standard_obj
    testthat::expect_equal(hash_computation_result,hash_gold_standard_obj$value)
  }else{
    hash_gold_standard_obj = list()
    hash_gold_standard_obj$value = hash_computation_result
    save(hash_gold_standard_obj,file = file_to_check)
  }
  
}
