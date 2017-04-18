#include <map>
#include <set>
#include <string>
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>

// Compile with: g++ -std=c++11 -o simulate_uc_events simulate_uc_events.cpp

// We represent sequences as strings of integers (from 0 to a possible 255).
// These represent one of a fixed number of chunks of the reference sequence.
// E.g. for DUP4 there are six breakpoints so 7 segments numbered 0 to 6.
// 0 - left of GYPE
// 1 - region around GYPE
// 2 - between GYPE and GYPB
// 3 - last exon of GYPB
// 4 - GYPB to left of GYPA
// 5 - left to mid GYPA
// 6 - GYPA onwards
//
// The reference chromosome is 0 1 2 3 4 5 6 while DUP4 is
// 0 1 2 1 5 4 5 4 5 6

typedef std::string ChunkSequence ;
struct RecombinationInfo ;
typedef std::map<
	ChunkSequence,
	std::pair<
		std::size_t,
		std::vector< RecombinationInfo >
	>
> CNVs ;
typedef CNVs::const_iterator Iterator ;

// This struct represents a recombination between two chromosomes.
// The recombinant copies left sequence up until just before pos1
// and then copies the right sequence from pos2 onwards.
// precondition: pos1 > 0 and pos1 < left.size(), pos2 > 0 and pos2 < right.size()
// (This corresponds to the requirement that the first and last chunks are never involved
// in deletion/duplication.)
struct RecombinationInfo {
	RecombinationInfo(
		Iterator const& left,
		std::size_t pos1,
		Iterator const& right,
		std::size_t pos2
	):
		m_left( left ),
		m_pos1( pos1 ),
		m_right( right ),
		m_pos2( pos2 )
	{
		assert( m_pos1 > 0 ) ;
		assert( m_pos2 > 0 ) ;
		assert( m_pos1 < m_left->first.size() ) ;
		assert( (m_pos2) < m_right->first.size() ) ;
	}

	RecombinationInfo( RecombinationInfo const& other ):
		m_left( other.m_left ),
		m_pos1( other.m_pos1 ),
		m_right( other.m_right ),
		m_pos2( other.m_pos2 )
	{}

	RecombinationInfo& operator=( RecombinationInfo const& other ) {
		m_left = other.m_left ;
		m_pos1 = other.m_pos1 ;
		m_right = other.m_right ;
		m_pos2 = other.m_pos2 ;
		return *this ;
	}

	// Generate the new recombinant sequence
	ChunkSequence generate() const {
		ChunkSequence result ;
		result.resize( m_pos1 + m_right->first.size() - m_pos2 ) ;
		std::copy( m_left->first.begin(), m_left->first.begin() + m_pos1, &result[0] ) ;
		std::copy( m_right->first.begin() + m_pos2, m_right->first.end(), &result[0] + m_pos1 ) ;
		return result ;
	}

	ChunkSequence const& left() const { return m_left->first ;}
	std::size_t const left_generation() const { return m_left->second.first ;}
	ChunkSequence const& right() const { return m_right->first ;}
	std::size_t const right_generation() const { return m_right->second.first ;}
	std::size_t const pos1() const { return m_pos1 ;}
	std::size_t const pos2() const { return m_pos2 ;}
private:
	Iterator m_left ;
	std::size_t m_pos1 ;
	Iterator m_right ;
	std::size_t m_pos2 ;
} ;

// Pretty-print a recombinationInfo.
// This prints the left and right recombining sides, with breakpoints.
// Also prints the generation of left and right sides.
std::ostream& operator<<( std::ostream& out, RecombinationInfo const& r ) {
	out << r.generate()
		<< "\t"
		<< (std::max( r.left_generation(), r.right_generation() ) + 1)
		<< "\t" ;
	for( std::size_t i = 0; i < r.left().size(); ++i ) {
		if( i == r.pos1() ) {
			out << "|" ;
		}
		out << r.left()[i] ;
	}
	out << "\t" ;
	for( std::size_t i = 0; i < r.right().size(); ++i ) {
		if( i == r.pos2() ) {
			out << "|" ;
		}
		out << r.right()[i] ;
	}
	
	out << "\t" << r.left()[r.pos1() - 1 ] << "|" << r.right()[r.pos2()] ;
	//out << "---" << "\"" << r.left() << "\":" << &r.left() << ":" << r.pos1() << ":" << "\"" << r.right() << "\":"<< &r.right()  << ":" << r.pos2() << ".\n" ;
	return out ;
}

// Return the list of histories that generate a given copy number profile
// (a profile = the counts of each chunk.  Rather than actually count, we
// represent each profile as a sorted sequence of chunks.)
void getHistoriesThatGenerateProfile(
	CNVs const& cnvs,
	ChunkSequence const& target,
	std::vector< RecombinationInfo >* result
) {
	for( Iterator i = cnvs.begin(); i != cnvs.end(); ++i ) {
		ChunkSequence sorted = i->first  ;
		std::sort( sorted.begin(), sorted.end() ) ;
		if( sorted == target ) {
			result->reserve( result->size() + i->second.second.size() ) ;
			result->insert( result->end(), i->second.second.begin(), i->second.second.end() ) ;
		}
	}
}

// Return the set of histories recorded as creating a given arrangement.
// (If the arrangement not present, return an empty list).
void getHistoriesThatGenerateArrangement(
	CNVs const& cnvs,
	ChunkSequence const& target,
	std::vector< RecombinationInfo >* result
) {
	Iterator i = cnvs.find( target ) ;
	if( i != cnvs.end() ) {
		result->insert( result->end(), i->second.second.begin(), i->second.second.end() ) ;
	}
}

// always return true
bool everything( std::size_t const, std::size_t const, ChunkSequence const& ) {
	return true ;
}

// always return true
bool earliest_generation( std::size_t const generation, ChunkSequence const&, std::size_t const cnv_generation ) {
	return generation == cnv_generation ;
}

// return true if the cnv matches the given profile
bool matchProfile( ChunkSequence cnv, ChunkSequence const& target ) {
	std::sort( cnv.begin(), cnv.end() ) ;
	return cnv == target ;
}

// return true if the cnv matches the given arrangement exactly
bool matchExactly( ChunkSequence cnv, ChunkSequence const& target ) {
	return cnv == target ;
}

template< typename Matcher >
CNVs generateCNVs(
	CNVs cnvs,
	std::size_t nGenerations,
	Matcher const& shouldRecordHistory = &everything
) {
	for( std::size_t generation = 1; generation <= nGenerations; ++generation ) {
		std::cerr << "Computing CNVS in generation " << generation << "...\n" ;
		CNVs new_cnvs ;
		std::size_t leftCount = 0 ;
		for( Iterator i = cnvs.begin(); i != cnvs.end(); ++i, ++leftCount ) {
			
			ChunkSequence const& left = i->first ;
			for( Iterator j = cnvs.begin(); j != cnvs.end(); ++j ) {
				ChunkSequence const& right = j->first ;
				for( std::size_t pos1 = 1; pos1 < left.size(); ++pos1 ) {
					for( std::size_t pos2 = 1; pos2 < right.size(); ++pos2 ) {
						RecombinationInfo r( i, pos1, j, pos2 ) ;
						ChunkSequence result = r.generate() ;
						CNVs::iterator where = new_cnvs.find( result ) ;
						if( where == new_cnvs.end() ) {
							where = new_cnvs.insert(
								std::make_pair(
									result,
									std::make_pair( generation, std::vector< RecombinationInfo >())
								)
							).first ;
						}
						if( shouldRecordHistory( generation, result, where->second.first ) ) {
							where->second.second.push_back( r ) ;
						}
						if( new_cnvs.size() % 1000000 == 0 ) {
							std::cerr << "Looked at " << leftCount << " of " << cnvs.size() << " CNVS on left.  " ;
							std::cerr << "(" << new_cnvs.size() << " CNVs and counting...)" << std::endl ;
						}
					}
				}
			}
		}
		cnvs.insert( new_cnvs.begin(), new_cnvs.end() ) ;
		std::cerr << "After generation " << generation << ": "
			<< cnvs.size() << " haplotypes.\n" ;
	}

	return cnvs ;
}

int main( int argc, char** argv ) {
	if( argc < 3 ) {
		std::cerr << "usage: simulate_uc_events <target> <number of generations>.\n" ;
		std::cerr
			<< "\nGiven a target sequence of chunks T and a number of generations n, this program\n"
			"generates all rearrangements made from a 'reference' haplotype (made up by linearly\n"
			"ordering chunks), by n generations of unequal crossover events.\n"
			"The first (leftmost in the target) and last (rightmost in the target) chunk\n"
			"are considered flanking chunks and never involved in recombination,\n"
			"other chunks must be between these in ASCII order.\n"
			"\nexample: simulate_uc_events 01215456 2.\n\n" ;
		return -1 ;
	}
	std::size_t const generations = std::atoi( argv[2] ) ;
	if( generations > 3 ) {
		std::cerr << "Oh dear: more than 3 generations may make your computer explode.\n" ;
		return -1 ;
	}
	ChunkSequence const target = argv[1] ;
	if( target.size() < 2 ) {
		std::cerr << "Oh dear, target should have at least 2 elements.\n" ;
		return -1 ;
	}
	if( target[0] >= target[target.size()-1] ) {
		std::cerr << "Oh dear, target first char should be less than last char (suggest using 0....n).\n" ;
		return -1 ;
	}
	
	// Create the target copy number profile
	// (we represent this by just sorting the target.)
	ChunkSequence sortedTarget = target ;
	std::sort( sortedTarget.begin(), sortedTarget.end() ) ;
	if( sortedTarget[0] != target[0] || sortedTarget[sortedTarget.size()-1] != target[target.size()-1]) {
		std::cerr << "Oh dear, target should contain only characters between the first and last (suggest using 0....n).\n" ;
	}

	// Construct the starting 'reference' sequence by taking all chunks in between the first and last chunk
	// in order.
	ChunkSequence startingSequence ;
	for( char c = target[0]; c <= target[target.size()-1]; ++c ) {
		startingSequence.push_back(c) ;
	}
	
	std::cout << "# simulate_uc_events\n"
		<< "# Starting sequence is: " << startingSequence << "\n"
		<< "# Target is: " << target << "\n"
		<< "# Attempting to make target in " << generations << " generations.\n" ;

	{
		CNVs start ;
		start[ startingSequence ] = std::make_pair( std::size_t(0), std::vector< RecombinationInfo >() );
		CNVs const& cnvs = generateCNVs(
			start,
			generations,
			// to save on memory, only record histories for cnvs matching the target profile
			// [&sortedTarget]( ChunkSequence const& a ) { return matchProfile( a, sortedTarget ) ; }
			//&everything
			//&earliest_generation
			[&sortedTarget,generations](
				std::size_t const generation,
				ChunkSequence const& cnv,
				std::size_t const cnv_generation
			) {
				return matchProfile( cnv, sortedTarget )
					|| ( generation == cnv_generation && generation < generations ) ;
			}
		) ;

		{
			std::vector< RecombinationInfo > histories ;
			getHistoriesThatGenerateProfile( cnvs, sortedTarget, &histories ) ;
			std::cout << "# After " << generations << " generations, a total of " << cnvs.size() << " CNVs are possible.\n" ;
			std::cout << "# After " << generations << " generations, a total of " << histories.size() << " histories match target copy number profile.\n" ;
			std::cout << "index\tresult\tgeneration\tleft\tright\tbreak\tcopy_number_profile_match\texact_match\n" ;

			std::set< std::string > components ;
			std::size_t count = 0 ;
			while( histories.size() > 0 ) {
				std::vector< RecombinationInfo > newHistories ;
				for( std::size_t i = 0; i < histories.size(); ++i ) {
					std::string result = histories[i].generate() ;
					bool exact_match = (result == target) ;
					std::sort( result.begin(), result.end() ) ;
					bool profile_match = (result == target) ;
					std::cout << count++ << "\t" << histories[i] << "\t"
						<< ( profile_match ? 1 : 0 ) << "\t"
						<< ( exact_match ? 1 : 0 ) << "\n" ;
					std::string const& left = histories[i].left() ;
					std::string const& right = histories[i].right() ;
					if( components.find( left ) == components.end() ) {
						getHistoriesThatGenerateArrangement( cnvs, left, &newHistories ) ;
						components.insert( left ) ;
					}
					if( components.find( right ) == components.end() ) {
						getHistoriesThatGenerateArrangement( cnvs, right, &newHistories ) ;
						components.insert( right) ;
					}
				}
				histories.swap( newHistories ) ;
			}
		}
	}
	std::cerr << "Thanks for using simulate_uc_events!\n" ;
}
