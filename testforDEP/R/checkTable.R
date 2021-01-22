

checkTable = function(Ts, n, testName){
  if(testName == "vexler"){
    cutoff_table = vexler_cutoff
  }

  else if(testName == "MIC"){
    cutoff_table = MIC_cutoff
  }

  else if(testName == "EL"){
    cutoff_table = EL_cutoff
  }


  low = floor(n/100)
  up = low + 1

  vec = ((1-((n-100*low)/100))*cutoff_table[low,] + ((n-100*low)/100)*cutoff_table[up,])
  return (sum(vec > Ts)/(length(vec)+1))

}
