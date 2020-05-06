expandGenCoordinate <- function(genomic.coordinate, env=parent.frame(), k) {
  # Expand coordinate to size k. If the coordinate is bigger than k, stop.
  #
  # genomic.coordinate    <string>      A variable name pointing to genomic coordinate <data.table>
  #                                     containing column "chromosome", "start", "end"
  # env                 <environment>   An environment object where the data.table exists.
  
  # Dependency
  #    Package: data.table
  
  if (class(genomic.coordinate)[1] != "character") {
    stop("Please input a variable name in string format.")
  }
  
  # check size
  if (env[[genomic.coordinate]][end-start+1 > k, .N] > 0) {
    message("The length of coordinates is more than k. The coordinates cannot be expanded.")
  } else if (env[[genomic.coordinate]][end-start+1 > k, .N] > 0) {
    message("Original coordinate is longer than DNA pattern length. Shrinking the coordinates...")
  }
  
  # shrink start and end coordinates to point to their midpoint
  env[[genomic.coordinate]][, `:=`(
    
    start = {
      
      len <- end-start+1
      even <- len %% 2 == 0
      odd <- len %% 2 == 1
      
      # if even length
      start[even] = start + len/2 - 1
      # if odd length
      start[odd] = start + len%/%2
      
      start
    },
    
    end = {
      # if even length
      end[even] = end - len/2 + 1
      # if odd length
      end[odd] = end - len%/%2
      end
    }
  )]
  
  # genomic coordinate length
  len <- env[[genomic.coordinate]][, unique(end - start + 1)]
  
  # if len is not either 1 or 2, something is wrong with the shrinking
  if (!len %in% 1:2 | length(len) > 1) {
    stop("Something is wrong. Contact the developer.")
  }
  
  # scale up to k; if middle point is odd number expand k/2 on both side; if even, expand k/2-1
  # if size k cannot make flanking size b/w the case pattern balance, change k to k+1.
  env[[genomic.coordinate]][, `:=`(
    start = {
      len <- end-start+1
      start[len == 2] <- start - as.integer(k/2+0.5) + 1
      start[len == 1] <- start - k%/%2
      start
    },
    end = {
      end[len == 2] <- end + as.integer(k/2+0.5) - 1
      end[len == 1] <- end + k%/%2
      end
    }
  )]
  
  if (env[[genomic.coordinate]][, unique(end-start+1) > k]) {
    
    message(paste0("--Uneven flank between the case region. Size of k is changed from ", k, " to ",
                   env[[genomic.coordinate]][, unique(end-start+1)], "."))
    message("--k is updated automatically.")
    env$k <- env[[genomic.coordinate]][, unique(end-start+1)]
    assign("k", env$k, envir = parent.frame())
  }
}