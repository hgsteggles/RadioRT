/*
 * ParseLua.hpp
 *
 *  Created on: 16 Jul 2015
 *      Author: harry
 */

#ifndef PARSELUA_HPP_
#define PARSELUA_HPP_

#include <string>

#include "Logger.hpp"

#include "lib/selene/include/selene.h"

std::string parseString(sel::Selector selector) {
	std::string str = selector;
	return str;
}

bool exists(sel::Selector selector) {
	std::string str = selector;

	return !(str.compare("") == 0);
}

template<typename T>
void parseLuaVariable(sel::Selector selector, T& var) {
	if (exists(selector)) {
		var = (T) selector;
	}
	else {
		Logger<ConsoleLogPolicy>::Instance().print<SeverityType::WARNING>("parseLuaVariable: <", selector.getName(), "> does not exist, using default.");
		Logger<FileLogPolicy>::Instance().print<SeverityType::WARNING>("parseLuaVariable: <", selector.getName(), "> does not exist, using default.");
	}
	Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>(selector.getName(), ": ", var);
}

template<>
void parseLuaVariable(sel::Selector selector, std::string& var) {
	if (exists(selector)) {
		var = parseString(selector);
	}
	else {
		Logger<ConsoleLogPolicy>::Instance().print<SeverityType::WARNING>("parseLuaVariable: <", selector.getName(), "> does not exist, using default.");
		Logger<FileLogPolicy>::Instance().print<SeverityType::WARNING>("parseLuaVariable: <", selector.getName(), "> does not exist, using default.");
	}
	Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>(selector.getName(), ": ", var);
}



#endif /* PARSELUA_HPP_ */
