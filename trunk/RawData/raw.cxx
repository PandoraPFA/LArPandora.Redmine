/// \file    raw.cxx
/// \brief   raw data utilities
/// \author  brebel@fnal.gov
/// \version $Id: raw.cxx,v 1.7 2010/04/16 13:52:25 brebel Exp $
#include "RawData/raw.h"
#include <iostream>
#include <string>
#include <bitset>
#include <cstdlib>

namespace raw {

  //----------------------------------------------------------
  void Compress(std::vector<short> &adc, raw::Compress_t compress)
  {
    if(compress == raw::kHuffman) CompressHuffman(adc);

    return;
  }

  //----------------------------------------------------------
  // if the compression type is kNone, copy the adc vector into the uncompressed vector
  void Uncompress(const std::vector<short> adc, std::vector<short> &uncompressed, raw::Compress_t compress)
  {
    if(compress == raw::kHuffman) UncompressHuffman(adc, uncompressed);
    else if(compress == raw::kNone){
      for(unsigned int i = 0; i < adc.size(); ++i) uncompressed[i] = adc[i];
    }

    return;
  }
  
  // the current Huffman Coding scheme used by uBooNE is
  // based on differences between adc values in adjacent time bins
  // the code is 
  // no change for 4 ticks --> 1
  // no change for 1 tick  --> 01
  // +1 change             --> 001
  // -1 change             --> 0001
  // +2 change             --> 00001
  // -2 change             --> 000001
  // +3 change             --> 0000001
  // -3 change             --> 00000001
  // abs(change) > 3       --> write actual value to short
  // use 15th bit to set whether a block is encoded or raw value
  // 1 --> Huffman coded, 0 --> raw
  // pad out the lowest bits in a word with 0's
  void CompressHuffman(std::vector<short> &adc)
  {
    std::vector<short> workvec(adc);
    std::vector<short> diffs(adc.size(), 0);
    diffs[0] = 0;

//     std::cout << adc.size() << " " << diffs.size() << std::endl;

    // loop over the rawadc and record the differences between 
    // successive entries.  the first entry is always left
    // as the raw adc value, just set the bit 15 to 0
    for(unsigned int i = 1; i < diffs.size(); ++i) diffs[i] = adc[i] - adc[i-1];

    // now loop over the diffs and do the Huffman encoding
    unsigned int cur = 1;
    unsigned int curb = 15;

    std::bitset<16> bset;
    bset.set(15);

    for(unsigned int i = 1; i < diffs.size(); ++i){
//       std::cout << i << " " << adc[i] << " current diff = " << diffs[i] << " " << curb << " " 
// 		<< bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 		<< std::endl;

      //check to see if the difference is too large that we have to put the entire adc value in
      if(abs(diffs[i]) > 3){

	// put the current value into the adc vec unless the current bit is 15, then there 
	// were multiple large difference values in a row
	if(curb != 15){
	  adc[cur] = bset.to_ulong();
	  ++cur;
	}

// 	std::cout << cur-1 << " " << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() << std::endl;
	bset.reset();
	bset.set(15);
	curb = 15;

	// put the current adcvalue in adc, with its bit 15 set to 0
// 	std::cout << "putting actual value  in " << workvec[i] << std::endl;
	if(workvec[i] > 0) adc[cur] = workvec[i];
	else{
	  std::bitset<16> tbit(-workvec[i]);
	  tbit.set(14);
// 	  std::cout << "negative value " << workvec[i] << " " 
// 		    << tbit.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 		    << std::endl;
	  adc[cur] = tbit.to_ulong();
	} 
	++cur;

      }
      else{
	// if the difference is 0, check to see what the next 3 differences are
	if(diffs[i] == 0){
	  if(i < diffs.size() - 3){
	    // if next 3 are also 0, set the next bit to be 1
	    if(diffs[i+1] == 0 && diffs[i+2] == 0 && diffs[i+3] == 0){
	      if(curb > 0){
		--curb;
		bset.set(curb);
// 		std::cout << "multiple zeros " << i << " " << i+3 << " " << diffs[i+1] << " " 
// 			  << " " << diffs[i+2] << " " << diffs[i+3] << " " << curb << " " 
// 			  << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 			  << std::endl;
		i += 3; 
	      }
	      else{	    
		adc[cur] = bset.to_ulong();
		++cur;
// 		std::cout << cur-1 << " " 
// 			  << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 			  << std::endl;
		// reset the bitset to be ready for the next word
		bset.reset();
		bset.set(15);
		bset.set(14); // account for the fact that this is a zero diff
		curb = 14;
		i += 3; 
	      } // end if curb is not big enough to put current difference in bset	  
	    } // end if next 3 are also zero
	    else{
	      // 0 diff is encoded as 01, so move the current bit one to the right
	      if(curb > 1){
		curb -= 2;
		bset.set(curb);
	      } // end if the current bit is large enough to set this one
	      else{	    
		adc[cur] = bset.to_ulong();
		++cur;
// 		std::cout << cur-1 << " " 
// 			  << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 			  << std::endl;
		// reset the bitset to be ready for the next word
		bset.reset();
		bset.set(15);
		bset.set(13); // account for the fact that this is a zero diff
		curb = 13;
	      } // end if curb is not big enough to put current difference in bset	  	    
	    } // end if next 3 are not also 0
	  }// end if able to check next 3
	  else{
	    // 0 diff is encoded as 01, so move the current bit one to the right
	    if(curb > 1){
	      curb -= 2;
	      bset.set(curb);
	    } // end if the current bit is large enough to set this one
	    else{	    
	      adc[cur] = bset.to_ulong();
	      ++cur;
// 	      std::cout << cur-1 << " " 
// 			<< bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 			<< std::endl;
	      // reset the bitset to be ready for the next word
	      bset.reset();
	      bset.set(15);
	      bset.set(13); // account for the fact that this is a zero diff
	      curb = 13;
	    } // end if curb is not big enough to put current difference in bset	  
	  }// end if not able to check the next 3
	}// end if current difference is zero
	else if( diffs[i] == 1){
	  if(curb > 2){
	    curb -= 3;
	    bset.set(curb);
	  }
	  else{
	    adc[cur] = bset.to_ulong();
	    ++cur;
// 	    std::cout << cur-1 << " " 
// 		      << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 		      << std::endl;
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(12); // account for the fact that this is a +1 diff
	    curb = 12;
	  } // end if curb is not big enough to put current difference in bset	  
	}// end if difference = 1
	else if( diffs[i] == -1){
	  if(curb > 3){
	    curb -= 4;
	    bset.set(curb);
	  }
	  else{
	    adc[cur] = bset.to_ulong();
	    ++cur;
// 	    std::cout << cur-1 << " " 
// 		      << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 		      << std::endl;
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(11); // account for the fact that this is a -1 diff
	    curb = 11;
	  } // end if curb is not big enough to put current difference in bset	  
	}// end if difference = -1
	else if( diffs[i] == 2){
	  if(curb > 4){
	    curb -= 5;
	    bset.set(curb);
	  }
	  else{
	    adc[cur] = bset.to_ulong();
	    ++cur;
// 	    std::cout << cur-1 << " " 
// 		      << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 		      << std::endl;
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(10); // account for the fact that this is a +2 diff
	    curb = 10;
	  } // end if curb is not big enough to put current difference in bset	  
	}// end if difference = 2
	else if( diffs[i] == -2){
	  if(curb > 5){
	    curb -= 6;
	    bset.set(curb);
	  }
	  else{
	    adc[cur] = bset.to_ulong();
	    ++cur;
// 	    std::cout << cur-1 << " " 
// 		      << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 		      << std::endl;
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(9); // account for the fact that this is a -2 diff
	    curb = 9;
	  } // end if curb is not big enough to put current difference in bset	  
	}// end if difference = -2
	else if( diffs[i] == 3){
	  if(curb > 6){
	    curb -= 7;
	    bset.set(curb);
	  }
	  else{
	    adc[cur] = bset.to_ulong();
	    ++cur;
// 	    std::cout << cur-1 << " " 
// 		      << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 		      << std::endl;
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(8); // account for the fact that this is a +3 diff
	    curb = 8;
	  } // end if curb is not big enough to put current difference in bset	  
	}// end if difference = 3
	else if( diffs[i] == -3){
	  if(curb > 7){
	    curb -= 8;
	    bset.set(curb);
	  }
	  else{
	    adc[cur] = bset.to_ulong();
	    ++cur;
// 	    std::cout << cur-1 << " " 
// 		      << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 		      << std::endl;
	    // reset the bitset to be ready for the next word
	    bset.reset();
	    bset.set(15);
	    bset.set(7); // account for the fact that this is a -3 diff
	    curb = 7;
	  } // end if curb is not big enough to put current difference in bset	  
	}// end if difference = -3
      }// end if difference is <= 3
    }// end loop over differences

//     std::cout << "ready to write last one " << cur << " " 
// 	      << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
// 	      << std::endl;

    //write out the last bitset
    adc[cur] = bset.to_ulong();
    ++cur;

    diffs.clear();
    workvec.clear();
 
    // resize the adc to be the number of entries used
    adc.resize(cur);

    return;
  }
  //--------------------------------------------------------
  // need to decrement the bit you are looking at to determine the deltas as that is how
  // the bits are set
  void UncompressHuffman(const std::vector<short> adc, std::vector<short> &uncompressed)
  {
//     std::cout << uncompressed.size() << std::endl;
    
    //the first entry in adc is a data value by construction
    uncompressed[0] = adc[0];

    unsigned int curu = 1;
    short curADC = uncompressed[0];

//     std::cout << "first ADC is " << curADC << std::endl;

    // loop over the entries in adc and uncompress them according to the
    // encoding scheme above the CompressHuffman method
    for(unsigned int i = 1; i < adc.size() && curu < uncompressed.size(); ++i){

      std::bitset<16> bset(adc[i]);

//       std::cout << curu << " " << adc.size() - i << " " << adc[i] << " " 
// 		<< bset.to_string< char,std::char_traits<char>,std::allocator<char> >() << std::endl;

      int numu = 0;

      //check the 15 bit to see if this entry is a full data value or not
      if( !bset.test(15) ){
	curADC = adc[i];
	if(bset.test(14)){
	  bset.set(14, false);
	  curADC = -1*bset.to_ulong();
	}
	uncompressed[curu] = curADC;

// 	std::cout << "bit 15 not set " << i << " " << adc[i] 
// 		  << " " << curADC << " " << curu << std::endl;

	++curu;
      }
      else{

//  	std::cout << "bit 15 set " << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() << std::endl;

	int  b       = 14;
	int  lowestb = 0;

	// ignore any padding with zeros in the lower order bits
	while( !bset.test(lowestb) && lowestb < 15) ++lowestb;

//  	std::cout << "zero padding skipped, lowest bit now at " << lowestb << std::endl;

	if(lowestb > 14){
	  std::cout << "encoded entry has no set bits!!! " << i << " "
		    << bset.to_string< char,std::char_traits<char>,std::allocator<char> >() 
		    << std::endl;
	  continue;
	}

	while( b >= lowestb){ 

	  // count the zeros between the current bit and the next on bit
	  int zerocnt = 0;
	  while( !bset.test(b-zerocnt) && b-zerocnt > lowestb) ++zerocnt;

	  b -= zerocnt;

// 	  std::cout << "bit now at " << b << " zerocnt = " << zerocnt 
// 		    << " lowest bit " << lowestb << " current uncompressed " << curu << std::endl;

	  if(zerocnt == 0){
	    for(int s = 0; s < 4; ++s){
	      uncompressed[curu] = curADC;
// 	      std::cout << curu << " " << uncompressed[curu] << std::endl;
	      ++curu;
	      ++numu;
	      if(curu > uncompressed.size()-1) break;
	    }
	    --b;
	  }
	  else if(zerocnt == 1){
	    uncompressed[curu] = curADC;
// 	    std::cout << curu << " " << uncompressed[curu] << std::endl;
	    ++curu;
	    ++numu;
	    --b;
	  }
	  else if(zerocnt == 2){
	    curADC += 1;
	    uncompressed[curu] = curADC;
// 	    std::cout << curu << " " << uncompressed[curu] << std::endl;
	    ++curu;
	    ++numu;
	    --b;
	  }
	  else if(zerocnt == 3){
	    curADC -= 1;
	    uncompressed[curu] = curADC;
// 	    std::cout << curu << " " << uncompressed[curu] << std::endl;
	    ++curu;
	    ++numu;
	    --b;
	  }
	  else if(zerocnt == 4){
	    curADC += 2;
	    uncompressed[curu] = curADC;
// 	    std::cout << curu << " " << uncompressed[curu] << std::endl;
	    ++curu;
	    ++numu;
	    --b;
	  }
	  else if(zerocnt == 5){
	    curADC -= 2;
	    uncompressed[curu] = curADC;
// 	    std::cout << curu << " " << uncompressed[curu] << std::endl;
	    ++curu;
	    ++numu;
	    --b;
	  }
	  else if(zerocnt == 6){
	    curADC += 3;
	    uncompressed[curu] = curADC;
// 	    std::cout << curu << " " << uncompressed[curu] << std::endl;
	    ++curu;
	    ++numu;
	    --b;
	  }
	  else if(zerocnt == 7){
	    curADC -= 3;
	    uncompressed[curu] = curADC;
// 	    std::cout << curu << " " << uncompressed[curu] << std::endl;
	    ++curu;
	    ++numu;
	    --b;
	  }

	  if(curu > uncompressed.size() - 1) break;

	}// end loop over bits
 
	if(curu > uncompressed.size() - 1) break;

//  	std::cout << "there were " << numu << " adc values uncompressed from this entry " << curu << std::endl;

      }// end if this entry in the vector is encoded

//       std::cout << curu << " " << i << " " << curADC << std::endl;

    }// end loop over entries in adc

    return;
  }

}
