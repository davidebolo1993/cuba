#ifndef PWALIGN_H
#define PWALIGN_H

#include <filesystem>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>


struct cmd_arguments_pwalign {
	std::vector<std::string> stringin{};
	std::string type{"global"};
	int match{2};
	int mismatch{-3};
	int gapopen{-4};
	int gapextend{-2};
};


void initialise_argument_parser_pwalign(seqan3::argument_parser & subparser, cmd_arguments_pwalign & args)
{
	subparser.info.description.push_back("Perform pairwise alignment between a couple of strings");
	subparser.add_positional_option(args.stringin, "a couple of strings");
	subparser.add_option(args.type, 'a', "alignment", "Alignment type. Choose between global or local.", seqan3::option_spec::DEFAULT, seqan3::value_list_validator{"global", "local"});
	subparser.add_option(args.match, 'm', "match", "Reward for a matching base", seqan3::option_spec::DEFAULT);
	subparser.add_option(args.mismatch, 'x', "mismatch", "Penalty for a mismatching base", seqan3::option_spec::DEFAULT);
	subparser.add_option(args.gapopen, 'g', "gapopen", "Penalty for opening a gap", seqan3::option_spec::DEFAULT);
	subparser.add_option(args.gapextend, 'e', "gapextend", "Penalty for extending a gap", seqan3::option_spec::DEFAULT);
};


int pwalign(seqan3::argument_parser & subparser)
{


	time_t my_time; 
	my_time= time(NULL);
	char *t = ctime(&my_time);
	cmd_arguments_pwalign args{};
	initialise_argument_parser_pwalign(subparser, args);

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

	if (args.stringin.size() != 2) {

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Error][" <<  t << "] A couple of strings must be provided. Provided " << args.stringin.size() << " sequences instead" << std::endl;
		return -1;

	}

	//initialize congifurations


	auto output_config = seqan3::align_cfg::output_score{} |
						 seqan3::align_cfg::output_begin_position{} |
						 seqan3::align_cfg::output_end_position{} |
						 seqan3::align_cfg::output_alignment{};

	auto config_global = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
														 seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
														 seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
														 seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}} |
				   seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{args.match}, seqan3::mismatch_score{args.mismatch}}} |
				   seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{args.gapopen}, seqan3::align_cfg::extension_score{args.gapextend}} |
				   output_config;


	auto config_local = seqan3::align_cfg::method_local{} |
					seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{args.match}, seqan3::mismatch_score{args.mismatch}}} |
					seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{args.gapopen}, seqan3::align_cfg::extension_score{args.gapextend}} |
					output_config;



	std::vector<seqan3::dna5> sequence1 {};
	std::vector<seqan3::dna5> sequence2 {};

	//convert strings to dna5

	for (char c : args.stringin.front()) sequence1.push_back(seqan3::assign_char_to(c, seqan3::dna5{})); //fill vector seq
	for (char c : args.stringin.back()) sequence2.push_back(seqan3::assign_char_to(c, seqan3::dna5{})); //fill vector seq

	//global alignment

	if (args.type == "global") {

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Performing global alignment" << std::endl;

		for (auto const & res : seqan3::align_pairwise(std::tie(sequence1, sequence2), config_global))

		{
			seqan3::debug_stream << "Alignment score: " << res.score() << std::endl;
			seqan3::debug_stream << "Sequence 1 alignment range: " << res.sequence1_begin_position()+1 << "," << res.sequence1_end_position() << std::endl;
			seqan3::debug_stream << "Sequence 2 alignment range: " << res.sequence2_begin_position()+1 << "," << res.sequence2_end_position() << std::endl;
			//seqan3::debug_stream << "Alignment:"  << std::endl << res.alignment() << std::endl;
			//std::cout << res.alignment()[0] << std::endl;
		}


	}

	else  { // is local

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Performing local alignment" << std::endl;

		for (auto const & res : seqan3::align_pairwise(std::tie(sequence1, sequence2), config_local))

		{
			seqan3::debug_stream << "Alignment score: " << res.score() << std::endl;
			seqan3::debug_stream << "Sequence 1 alignment range: " << res.sequence1_begin_position()+1 << "," << res.sequence1_end_position() << std::endl;
			seqan3::debug_stream << "Sequence 2 alignment range: " << res.sequence2_begin_position()+1 << "," << res.sequence2_end_position() << std::endl;
			//seqan3::debug_stream << "Alignment:"  << std::endl << res.alignment() << std::endl;
		}


	}


	t = ctime(&my_time);
	t[strlen(t)-1] = '\0';
	std::cout << "[Message][" <<  t << "] Done" << std::endl;

	return 0;

}

#endif
