#ifndef JSON_LOGGER_H
#define JSON_LOGGER_H

#include <iostream>
#include <string>

#include "constants.h"

enum EVENT_TYPE {
    PRELUDE = -1,
    SOURCE,
    DESTINATION,
    GENERATING,
    EXPANDING,
    CLOSING,
    UPDATING,
    END,
    FINAL

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

        void log(enum EVENT_TYPE event_type, warthog::sn_id_t node_id = -1){
            int32_t x, y;
            warthog::search_node* node;
            switch(event_type){
                case PRELUDE:
                    std::cout << "{\"nodeStructure\": [{\"type\": \"rectangle\", \"variables\": {\"x\":"
                        << " \"x\", \"y\": \"y\"}, \"persisted\": true, \"drawPath\": true}]," << std::endl;
                    std::cout << "\"eventList\":[" << std::endl;
                    break;
                case SOURCE:
                case DESTINATION:
                    expander_->get_xy(node_id, x, y);
                    std::cout << "{\"type\":\"" << event_strings[event_type] << "\",\"id\":" << node_id << ",\"variables\":{\"x\":"
                        << x << ",\"y\":" << y << "}}," << std::endl;
                    break;
                case EXPANDING:
                case GENERATING:
                case CLOSING:
                case UPDATING:
                case END:
                    expander_->get_xy(node_id, x, y);
                    node = expander_->generate(node_id);
                    std::cout << "{\"type\":\"" << event_strings[event_type] << "\",\"id\":" << node_id << ",\"variables\":{\"x\":"
                        << x << ",\"y\":" << y << "},\"g\":" << node->get_g() << ",\"f\":" << node->get_f()
                        << ",\"pId\":" << node->get_parent() << "}";
                    if(event_type == END){
                        std::cout << "\n]}" << std::endl;
                    } else {
                        std::cout << "," << std::endl;
                    }
                    break;
                case FINAL:
                    //final path logging can be implemented here if required.
                    break;
                default:
                    std::cout << "UNHANDLED" << std::endl;
                    break;

            }
        }
    
    private:
        E* expander_;
        std::string event_strings[7] = {
        "source",
        "destination",
        "generating",
        "expanding",
        "closing",
        "updating",
        "end"
    };
};

}

#endif