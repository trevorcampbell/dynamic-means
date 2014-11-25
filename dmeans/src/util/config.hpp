#ifndef __CONFIG_HPP
#include<map>
#include<string>
#include<iostream>
#include<sstream>

class Config{
	public:
		enum Type{
			REQUIRED,
			OPTIONAL
		};
		template<class T>
			void set(std::string key, T val){
				std::ostringstream oss;
				oss << val;
				if (!oss.good()){
					throw ToStringConversionException(key);
				}
				params[key] = oss.str();
			}
		template<class T>
			T get(std::string key, Type t, T default_value){
				if(params.find(key) != params.end()){
					T res;
					std::istringstream iss;
					iss.str(params[key]);
					iss >> res;
					if (!iss.good()){
						throw FromStringConversionException(key, params[key]);
					}
					return res;
				} else if (t == Type::OPTIONAL){
					return default_value;
				} else {
					throw ParameterNotFoundException(key);
				}
			}
	private:
		std::map<std::string, std::string> params;

		class ParameterNotFoundException{
			ParameterNotFoundException(std::string s){
				std::cout << "Could not find parameter " << s << " in the config object." << std::endl;
			}
		};
		class FromStringConversionException{
			FromStringConversionException(std::string k, std::string s){
				std::cout << "Could not convert parameter " << k << " = " << s << " from string to a usable type for output." << std::endl;
			}
		};
		class ToStringConversionException{
			ToStringConversionException(std::string k){
				std::cout << "Could not store parameter " << k << " in the config -- the value could not be converted to string." << std::endl;
			}
		};
};

#define __CONFIG_HPP
#endif /* __CONFIG_HPP */
