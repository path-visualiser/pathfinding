#ifndef JSON_LOGGER_H
#define JSON_LOGGER_H

#include <iostream>
#include <string>
#include <array>

#include "search_node.h"

enum EVENT_TYPE {
    PRELUDE = -1,
    SOURCE,
    DESTINATION,
    GENERATING,
    EXPANDING,
    CLOSING,
    UPDATING,
    END

};



namespace warthog
{

// E is the expansion policy
template <class E>
class json_logger {
    public:
        json_logger(E* expander){
            expander_ = expander;
        }

        void log(enum EVENT_TYPE event_type, warthog::search_node* node){
            std::cout << expander_->log_string(event_type, node);
        }

        void final_path(){
            std::cout << "{\"paths\": []}" << std::endl
                << "]}" << std::endl;
        }
    
    private:
        E* expander_;


};

}

#endif