#ifndef FIND_H
#define FIND_H

#include <filesystem>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/fm_index.hpp> //for using the fm_index
#include <seqan3/search/fm_index/bi_fm_index.hpp>  //for using the bi_fm_index
#include <cereal/archives/binary.hpp>
#include <seqan3/search/search.hpp> //for searching
#include <seqan3/alphabet/all.hpp>

struct cmd_arguments_find {
	std::string stringin;
	std::string filein;
	int maxerr {0};
	bool bidirectional {true};
	bool all {true};
};


void initialise_argument_parser_find(seqan3::argument_parser & subparser, cmd_arguments_find & args)
{
	subparser.info.description.push_back("Search for a string in a (bidirectional) fm-index");
	subparser.add_positional_option(args.stringin, "input string to search for"); 
	subparser.add_flag(args.bidirectional, 'b', "bidirectional", "the index is bidirectional (is a .bifmi file)",seqan3::option_spec::DEFAULT);
	subparser.add_option(args.filein, 'f', "fmindex", "input (bi-)fm-index", seqan3::option_spec::REQUIRED);
	subparser.add_option(args.maxerr, 'e', "error", "maximum number of errors for approximate search", seqan3::option_spec::DEFAULT);
	subparser.add_flag(args.all, 'a', "all", "report all the hits",seqan3::option_spec::DEFAULT);
};


void bi_fmi_matcher(seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::collection> & indexin, std::vector<seqan3::dna5> & query, auto & cfg)
{

	//seqan3::debug_stream << search(query, indexin,cfg) << std::endl;
	

	auto search_results = search(query, indexin, cfg);

	//std::cout << "Found " << search_results.size() << " matches" >> std::endl;

	int counter=0;

	for (auto && hit : search_results) {

		seqan3::debug_stream << "Hit found on sequence " << hit.reference_id() +1 << ", starting at base " << hit.reference_begin_position() +1  << std::endl;
		counter ++;
	}

	if (counter == 0) {

		std::cout << "No hit found" << std::endl;

	}

};




void fmi_matcher(seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> & indexin, std::vector<seqan3::dna5> & query, auto & cfg)
{

	auto search_results = search(query, indexin, cfg);

	//std::cout << "Found " << search_results.size() << " matches" >> std::endl;

	int counter=0;

	for (auto && hit : search_results) {

		seqan3::debug_stream << "Hit found on sequence " << hit.reference_id() +1 << ", starting at base " << hit.reference_begin_position() +1  << std::endl;
		counter ++;
	}

	if (counter == 0) {

		std::cout << "No hit found" << std::endl;

	}

};



int find(seqan3::argument_parser & subparser)
{


	time_t my_time; 
	my_time= time(NULL);
	char *t = ctime(&my_time);
	cmd_arguments_find args{};
	initialise_argument_parser_find(subparser, args);

	try
	{
		subparser.parse();
	}
	
	catch (seqan3::argument_parser_error const & ext)
	{
		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Error][" <<  t << "] Wrong command-line argument" << std::endl;
		seqan3::debug_stream << ext.what() << std::endl;
		return -1;
	}

	std::vector<seqan3::dna5> sequence {};
	std::string fin;

	seqan3::search_cfg::hit hit_dynamic{seqan3::search_cfg::hit_all_best{}};
	if (args.all) hit_dynamic = seqan3::search_cfg::hit_all{};

	seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{args.maxerr}} | hit_dynamic;

	if (args.bidirectional) {

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Loading bidirectional fm-index" << std::endl;
		seqan3::bi_fm_index<seqan3::dna5, seqan3::text_layout::collection> indexin;
		fin=std::filesystem::canonical(args.filein).string();

		if (fin.substr(fin.find_last_of(".") + 1) != "bifmi") {

			t = ctime(&my_time);
			t[strlen(t)-1] = '\0';
			std::cout << "[Error][" <<  t << "] Wrong filename extension to bi-fm-index. If this is a .fmi file, remove the -b/--bidirectional flag" << std::endl;
			return -1;
		
		} // extension is wrong, stop

		{
		std::ifstream is{fin, std::ios::binary};
		cereal::BinaryInputArchive iarchive{is};
		iarchive(indexin);
		}

		for (char c : args.stringin) sequence.push_back(seqan3::assign_char_to(c, seqan3::dna5{})); //fill vector seq

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Searching through the bidirectional fm-index" << std::endl;
		bi_fmi_matcher(indexin, sequence,cfg);


	
	} else {


		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Loading fm-index" << std::endl;
		seqan3::fm_index<seqan3::dna5, seqan3::text_layout::collection> indexin;
		fin=std::filesystem::canonical(args.filein).string();

		if (fin.substr(fin.find_last_of(".") + 1) != "fmi") {

			t = ctime(&my_time);
			t[strlen(t)-1] = '\0';
			std::cout << "[Error][" <<  t << "] Wrong filename extension to fm-index. If this is a .bifmi file, add the -b/--bidirectional flag" << std::endl;
			return -1;
		
		} // extension is wrong, stop

		{
		std::ifstream is{fin, std::ios::binary};
		cereal::BinaryInputArchive iarchive{is};
		iarchive(indexin);
		}

		for (char c : args.stringin) sequence.push_back(seqan3::assign_char_to(c, seqan3::dna5{})); //fill vector seq
		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Searching through the fm-index" << std::endl;
		fmi_matcher(indexin, sequence, cfg);

	}


	t = ctime(&my_time);
	t[strlen(t)-1] = '\0';
	std::cout << "[Message][" <<  t << "] Done" << std::endl;

	return 0;

}

#endif
