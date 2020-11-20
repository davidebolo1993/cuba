#ifndef INDEX_H
#define INDEX_H

#include <zlib.h>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <fstream>
#include <filesystem>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/fm_index.hpp> // for using the fm_index
#include <seqan3/search/fm_index/bi_fm_index.hpp>  // for using the bi_fm_index
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp> 
#include <seqan3/alphabet/all.hpp>

//headers
#include "kseq.h"

struct cmd_arguments_index {
	std::vector<std::string> filein{};
	std::string fileout {"out.fmi"};
	std::string vecout;
	bool bidirectional {true};
};


KSEQ_INIT(gzFile, gzread)

void initialise_argument_parser_index(seqan3::argument_parser & subparser, cmd_arguments_index & args)
{
	subparser.info.description.push_back("Create a full-searchable (bidirectional) fm-index from fasta/fastq file/s");
	subparser.add_positional_option(args.filein, "input fastq/fasta file/s, optionally gzip-compressed"); 
	subparser.add_flag(args.bidirectional, 'b', "bidirectional", "create bidirectional fm-index (out.bifmi)",seqan3::option_spec::DEFAULT);
	subparser.add_option(args.fileout, 'f', "fmindex", "output (bidirectional) fm-index", seqan3::option_spec::DEFAULT);
	subparser.add_option(args.vecout, 'v', "vector", "serialize vector of sequences to file");
};


int index(seqan3::argument_parser & subparser)
{

	time_t my_time; 
	my_time= time(NULL);
	char *t = ctime(&my_time);
	cmd_arguments_index args{};
	initialise_argument_parser_index(subparser, args);

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

	gzFile fp;
	kseq_t *seq;
	int l;
	std::vector<seqan3::dna5_vector> sequences = {};
	std::string fmout = std::filesystem::absolute(std::filesystem::weakly_canonical(args.fileout).string()).string();
	std::string tmpfile;
	if (!args.vecout.empty()) tmpfile = std::filesystem::absolute(std::filesystem::weakly_canonical(args.vecout).string()).string();

	for (std::vector<std::string>::const_iterator i = args.filein.begin(); i != args.filein.end(); ++i) {

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Reading " << std::filesystem::canonical(*i) << std::endl;

		fp = gzopen(std::filesystem::canonical(*i).string().c_str(), "r");
		seq = kseq_init(fp);
		std::string n,s,o,q;

		while ((l = kseq_read(seq)) >= 0) {

			std::vector<seqan3::dna5> sequence {};
			n=seq->name.s;
			s=seq->seq.s;
			
			//if (seq->qual.l) { //is fastq

				//o=seq->comment.s;
				//q=seq->qual.s;
			//} //this is not used for the time being, but we can add here filters based on quality

			for (char c : s) sequence.push_back(seqan3::assign_char_to(c, seqan3::dna5{})); //fill vector seq
			sequences.push_back(sequence);
		}

		gzclose(fp);

	}

	kseq_destroy(seq);

	if (!tmpfile.empty()) {

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Storing sequences vector to file" << std::endl;

		{
		std::ofstream vecout(tmpfile, std::ios::binary);
		cereal::BinaryOutputArchive archive(vecout);
		archive(sequences); 
		}

	}

	if (args.bidirectional) {

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Building bidirectional fm-index" << std::endl;
		seqan3::bi_fm_index indexout{sequences};

		if (fmout.substr(fmout.find_last_of(".") + 1) != "bifmi") {

			t = ctime(&my_time);
			t[strlen(t)-1] = '\0';
			std::cout << "[Warning][" <<  t << "] Wrong filename extension to bi-fm-index. Changing to .bifmi" << std::endl;
			fmout.replace(fmout.find_last_of(".") + 1, fmout.length()-fmout.find_last_of(".")+1, "bifmi");
		
		} // extension is wrong, replace

		{
		std::ofstream os{fmout, std::ios::binary};
		cereal::BinaryOutputArchive oarchive{os};
		oarchive(indexout);
		}

	} else {

		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Message][" <<  t << "] Building fm index" << std::endl;
		seqan3::fm_index indexout{sequences};

		if (fmout.substr(fmout.find_last_of(".") + 1) != "fmi") {

			t = ctime(&my_time);
			t[strlen(t)-1] = '\0';
			std::cout << "[Warning][" <<  t << "] Wrong filename extension to fm-index. Changing to .fmi" << std::endl;
			fmout.replace(fmout.find_last_of(".") + 1, fmout.length()-fmout.find_last_of(".")+1, "fmi");

		} // extension is wrong, replace

		{
		std::ofstream os{fmout, std::ios::binary};
		cereal::BinaryOutputArchive oarchive{os};
		oarchive(indexout);
		}
	}

	t = ctime(&my_time);
	t[strlen(t)-1] = '\0';
	std::cout << "[Message][" <<  t << "] Done" << std::endl;


	return 0;
}

#endif
