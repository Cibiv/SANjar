LS1 <- data.table(
  day  =c(  1,  1,  1,  2,  2,  2,  2,  2,   3),
  sid  =c('a','a','b','a','a','b','b','c', 'd'),
  reads=c(  2,  3,  4,  5,  6,  7,  8,  9,  10)
)

OS1 <- data.table(
  day  =c(  1,  1,  2,  2,  2),
  cells=c( 10, 12, 22, 23, 24)
)

test_that('LTData', {
  lt <- LTData(lineagesizes=LS1, organoidsizes=OS1)
  expect_equal(lt$sequencing$day, c(1,1,2,2,2,3))
  expect_equal(lt$sequencing$sid, c('a','b','a','b','c','d'))
  
  lt <- estimate_sequencing_parameters(lt, pcr_efficiency=0.5)
  expect_equal(lt$sequencing$day, c(1,1,2,2,2,3))
  expect_equal(lt$sequencing$sid, c('a','b','a','b','c','d'))
  expect_equal(lt$sequencing$library_size, c(5,4,11,15,9,10))
  expect_equal(lt$sequencing$pcr_efficiency, rep(0.5, 6))
  expect_equal(lt$sequencing$phantom_threshold, c(2,4,5,7,9,10))
  
  lt <- normalize_library_sizes(lt)
  expect_equal(lt$sequencing$day, c(1,1,2,2,2,3))
  expect_equal(lt$sequencing$sid, c('a','b','a','b','c','d'))
  expect_equal(lt$sequencing$library_size, c(4,4,9,9,9,10))
  expect_equal(lt$sequencing$pcr_efficiency, rep(0.5, 6))
  expect_equal(lt$sequencing$phantom_threshold, c(4,4,9,9,9,10))
  expect_equal(lt$lineagesizes$day, c(1,1,1,2,2,2,2,2,3))
  expect_equal(lt$lineagesizes$sid, c('a','a','b','a','a','b','b','c', 'd'))
  # Scale reads according to library sizes and re-apply threshold 
  r <- c(2,3,4,5,6,7,8,9,10) * (rep(c(4,4,9,9,9,10), times=c(2,1,2,2,1,1)) / 
                                rep(c(5,4,11,15,9,10), times=c(2,1,2,2,1,1)))
  r <- ifelse(r >= rep(c(4,4,9,9,9,10),times=c(2,1,2,2,1,1)), r, 0)
  expect_equal(lt$lineagesizes$reads, r)
})

test_that('LTData with groups', {
  expect_error(LTData(lineagesizes=rbind(cbind(g=1, LS1),
                                         cbind(g=2, LS1)),
                      organoidsizes=rbind(cbind(g=1, OS1)),
                      groups="g"),
               "organoidsizes table is missing an entry for g=2")
  
  lt <- LTData(lineagesizes=rbind(cbind(g=1, LS1),
                                   cbind(g=2, LS1[day<=2, list(day, sid, reads=2*reads)])),
                organoidsizes=rbind(cbind(g=1, OS1),
                                    cbind(g=2, OS1)),
                groups="g")
  expect_equal(lt$sequencing[g==1, day], c(1,1,2,2,2,3))
  expect_equal(lt$sequencing[g==1, sid], c('a','b','a','b','c','d'))
  expect_equal(lt$sequencing[g==2, day], c(1,1,2,2,2))
  expect_equal(lt$sequencing[g==2, sid], c('a','b','a','b','c'))
  
  lt <- estimate_sequencing_parameters(lt, pcr_efficiency=0.5)
  expect_equal(lt$sequencing[g==1, day], c(1,1,2,2,2,3))
  expect_equal(lt$sequencing[g==1, sid], c('a','b','a','b','c','d'))
  expect_equal(lt$sequencing[g==1, library_size], c(5,4,11,15,9,10))
  expect_equal(lt$sequencing[g==1, pcr_efficiency], rep(0.5, 6))
  expect_equal(lt$sequencing[g==1, phantom_threshold], c(2,4,5,7,9,10))
  expect_equal(lt$sequencing[g==2, day], c(1,1,2,2,2))
  expect_equal(lt$sequencing[g==2, sid], c('a','b','a','b','c'))
  expect_equal(lt$sequencing[g==2, library_size], c(10,8,22,30,18))
  expect_equal(lt$sequencing[g==2, pcr_efficiency], rep(0.5, 5))
  expect_equal(lt$sequencing[g==2, phantom_threshold], c(4,8,10,14,18))
  expect_equal(lt$sequencing[g==2, day], c(1,1,2,2,2))
  expect_equal(lt$sequencing[g==2, sid], c('a','b','a','b','c'))
  
  lt <- normalize_library_sizes(lt)
  expect_equal(lt$sequencing[g==1, day], c(1,1,2,2,2,3))
  expect_equal(lt$sequencing[g==1, sid], c('a','b','a','b','c','d'))
  expect_equal(lt$sequencing[g==1, library_size], c(4,4,9,9,9,10))
  expect_equal(lt$sequencing[g==1, pcr_efficiency], rep(0.5, 6))
  expect_equal(lt$sequencing[g==1, phantom_threshold], c(4,4,9,9,9,10))
  expect_equal(lt$lineagesizes[g==1, day], c(1,1,1,2,2,2,2,2,3))
  expect_equal(lt$lineagesizes[g==1, sid], c('a','a','b','a','a','b','b','c', 'd'))
  # Scale reads according to library sizes and re-apply threshold 
  r <- c(2,3,4,5,6,7,8,9,10) * (rep(c(4,4,9,9,9,10), times=c(2,1,2,2,1,1)) / 
                                  rep(c(5,4,11,15,9,10), times=c(2,1,2,2,1,1)))
  r <- ifelse(r >= rep(c(4,4,9,9,9,10),times=c(2,1,2,2,1,1)), r, 0)
  expect_equal(lt$lineagesizes[g==1, reads], r)
  expect_equal(lt$sequencing[g==2, day], c(1,1,2,2,2))
  expect_equal(lt$sequencing[g==2, sid], c('a','b','a','b','c'))
  expect_equal(lt$sequencing[g==2, library_size], c(8,8,18,18,18))
  expect_equal(lt$sequencing[g==2, pcr_efficiency], rep(0.5, 5))
  expect_equal(lt$sequencing[g==2, phantom_threshold], c(8,8,18,18,18))
  expect_equal(lt$lineagesizes[g==2, day], c(1,1,1,2,2,2,2,2))
  expect_equal(lt$lineagesizes[g==2, sid], c('a','a','b','a','a','b','b','c'))
  # Scale reads according to library sizes and re-apply threshold 
  r <- c(4,6,8,10,12,14,16,18) * (rep(c(8,8,18,18,18), times=c(2,1,2,2,1)) / 
                                  rep(c(10,8,22,30,18), times=c(2,1,2,2,1)))
  r <- ifelse(r >= rep(c(8,8,18,18,18),times=c(2,1,2,2,1)), r, 0)
  expect_equal(lt$lineagesizes[g==2, reads], r)
})
