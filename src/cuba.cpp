#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 

//headers
#include "index.h"
#include "find.h"
#include "pwalign.h"


inline void asciiArt() {

	std::cout << std::endl;
	std::cout << std::endl;
	std::cout <<" ________  ___  ___  ________  ________     "<< std::endl;
	std::cout <<"|\\   ____\\|\\  \\|\\  \\|\\   __  \\|\\   __  \\    "<< std::endl;
	std::cout <<"\\ \\  \\___|\\ \\  \\\\\\  \\ \\  \\|\\ /\\ \\  \\|\\  \\   "<< std::endl;
	std::cout <<" \\ \\  \\    \\ \\  \\\\\\  \\ \\   __  \\ \\   __  \\  "<< std::endl;
	std::cout <<"  \\ \\  \\____\\ \\  \\\\\\  \\ \\  \\|\\  \\ \\  \\ \\  \\ "<< std::endl;
	std::cout <<"   \\ \\_______\\ \\_______\\ \\_______\\ \\__\\ \\__\\"<< std::endl;
	std::cout <<"    \\|_______|\\|_______|\\|_______|\\|__|\\|__|"<< std::endl;                                            
	std::cout << std::endl;
	std::cout << std::endl;
}


int main(int argc, char const ** argv)
{
	seqan3::argument_parser top_level_parser{"cuba",
											 argc,
											 argv,
											 seqan3::update_notifications::off,
											 {"index", "find", "pwalign"}};

	// Top level parser
	top_level_parser.info.description.push_back("A collection of C++ modules based on ... to handle ... data efficiently");
	top_level_parser.info.short_description="cuba: C++ Utilities for BioinformAtics";
	top_level_parser.info.author = "Davide Bolognini";
	top_level_parser.info.version = "1.0";
	top_level_parser.info.email = "davidebolognini7@gmail.com";
	top_level_parser.info.url = "https://github.com/davidebolo1993/cuba";
	//parse top level parser
	time_t my_time; 
	my_time= time(NULL);
	char *t = ctime(&my_time);
	asciiArt();

	try
	{
		top_level_parser.parse(); // trigger command line parsing
	}

	catch (seqan3::argument_parser_error const & ext) // catch errors
	{
		t = ctime(&my_time);
		t[strlen(t)-1] = '\0';
		std::cout << "[Error][" <<  t << "] Wrong command-line argument" << std::endl;
		seqan3::debug_stream << ext.what() << std::endl; //check error
		return -1;
	}

	seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser

	t = ctime(&my_time);
	t[strlen(t)-1] = '\0';

	if (sub_parser.info.app_name == std::string_view{"cuba-index"}) {
		std::cout << "[Message][" <<  t << "] cuba index" << std::endl;
		return index(sub_parser);
	} 
	else if ( sub_parser.info.app_name == std::string_view{"cuba-find"}) {
		std::cout << "[Message][" <<  t << "] cuba find" << std::endl;
		return find(sub_parser);
	}
	else if ( sub_parser.info.app_name == std::string_view{"cuba-pwalign"}) {
		std::cout << "[Message][" <<  t << "] cuba align" << std::endl;
		return pwalign(sub_parser);
	}	
	return 0;
}
