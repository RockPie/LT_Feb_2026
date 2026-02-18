#include "H2GCROC_Common.hxx"
#include <iostream>
#include <unordered_set>
#include <unistd.h>
#include "TCanvas.h" 
#include "TVectorD.h"
#include "TVector.h"
#include "TF1.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TApplication.h"
#include "TGaxis.h"
#include "easylogging++.h"
#include "argparse/argparse.hpp"

#define FPGA_CHANNEL_NUMBER 152
#define HEARTBEAT_HEADER_START 0x23
#define HEARTBEAT_HEADER_END   0x23
#define DATA_PACKET_HEADER_START 0x23
#define DATA_PACKET_HEADER_END   0x23
#define LINE_NUMBER_PER_FPGA  20

INITIALIZE_EASYLOGGINGPP

int main(int argc, char **argv){

    // * --- Initial settings -----------------------------------------------------------
    // * --------------------------------------------------------------------------------
    std::string script_input_file, script_output_file;
    int script_n_events;
    std::string script_name = __FILE__;
    std::string script_version = "0.1";
    std::string script_output_folder;

    START_EASYLOGGINGPP(argc, argv);
    set_easylogger();

    script_name = script_name.substr(script_name.find_last_of("/\\") + 1).substr(0, script_name.find_last_of("."));

    argparse::ArgumentParser program(script_name, script_version);

    program.add_argument("-f", "--file").help("Input .h2g file").required();
    program.add_argument("-o", "--output").help("Output .root file").required();
    program.add_argument("-e", "--events").help("Number of events to process").default_value(std::string("-1"));
    try {
        program.parse_args(argc, argv);
        script_input_file  = program.get<std::string>("--file");
        script_output_file = program.get<std::string>("--output");
        auto script_n_events_str = program.get<std::string>("--events");
        script_n_events    = std::stoi(script_n_events_str);
    } catch (const std::runtime_error& err) {
        LOG(ERROR) << err.what();
        LOG(INFO) << program;
        return 1;
    }

    LOG(INFO) << "Input file: " << script_input_file;

    if (access(script_input_file.c_str(), F_OK) == -1) {
        LOG(ERROR) << "Input file " << script_input_file << " does not exist!";
        return 1;
    }
    if (script_output_file.substr(script_output_file.find_last_of(".") + 1) != "root") {
        LOG(ERROR) << "Output file " << script_output_file << " should end with .root!";
        return 1;
    }
    script_output_folder = script_output_file.substr(0, script_output_file.find_last_of("/\\"));
    if (script_output_folder.empty()) {
        script_output_folder = "./dump/" + script_name;
    }
    if (access(script_output_folder.c_str(), F_OK) == -1) {
        LOG(INFO) << "Creating output folder " << script_output_folder;
        if (mkdir(script_output_folder.c_str(), 0777) == -1) {
            LOG(ERROR) << "Failed to create output folder " << script_output_folder;
            return 1;
        }
    }
    if (access(script_output_file.c_str(), F_OK) != -1) {
        LOG(WARNING) << "Output file " << script_output_file << " already exists!";
    }

    LOG(INFO) << "Script name: " << script_name;
    LOG(INFO) << "Input file: " << script_input_file;
    LOG(INFO) << "Output file: " << script_output_file << " in " << script_output_folder;
    LOG(INFO) << "Number of events: " << script_n_events;

    // * --- Create output file ---------------------------------------------------------
    // * --------------------------------------------------------------------------------
    TFile *output_root = new TFile(script_output_file.c_str(), "RECREATE");
    if (output_root->IsZombie()) {
        LOG(ERROR) << "Failed to create output file " << script_output_file;
        return 1;
    }
    TTree* output_tree = new TTree("data_tree", "Data tree");
    output_tree->SetDirectory(output_root);

    // create branches
    // fpga_id is a 16-bit unsigned integer
    // timestamp is a 64-bit unsigned integer
    // daqh_list is 4 32-bit unsigned integers
    // tc and tp are boolen values
    auto *branch_fpga_id    = new UShort_t();   // 0 - 8
    auto *branch_timestamp  = new ULong64_t();  // 64 bits
    auto *branch_daqh_list  = new UInt_t[4];    // 32 bits
    auto *branch_tc_list    = new Bool_t[FPGA_CHANNEL_NUMBER];
    auto *branch_tp_list    = new Bool_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val0_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val1_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_val2_list  = new UInt_t[FPGA_CHANNEL_NUMBER];
    auto *branch_crc32_list = new UInt_t[4];
    auto *branch_last_heartbeat = new UInt_t();

    output_tree->Branch("fpga_id", branch_fpga_id, "fpga_id/s");
    output_tree->Branch("timestamp", branch_timestamp, "timestamp/l");
    output_tree->Branch("daqh_list", branch_daqh_list, "daqh_list[4]/i");
    output_tree->Branch("tc_list", branch_tc_list, ("tc_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/O").c_str());
    output_tree->Branch("tp_list", branch_tp_list, ("tp_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/O").c_str());
    output_tree->Branch("val0_list", branch_val0_list, ("val0_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    output_tree->Branch("val1_list", branch_val1_list, ("val1_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    output_tree->Branch("val2_list", branch_val2_list, ("val2_list[" + std::to_string(FPGA_CHANNEL_NUMBER) + "]/i").c_str());
    output_tree->Branch("crc32_list", branch_crc32_list, "crc32_list[4]/i");
    output_tree->Branch("last_heartbeat", branch_last_heartbeat, "last_heartbeat/i");

    // * --- Read input file ------------------------------------------------------------
    // * --------------------------------------------------------------------------------
    int legal_line_id_list[] = {0x00, 0x01, 0x02, 0x03, 0x04};
    int legal_asic_id_list[] = {0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7};
    int legal_half_id_list[] = {0x24, 0x25};
    std::vector <int> legal_fpga_id_list;
    legal_fpga_id_list.reserve(32);

    const int input_file_read_chunk_size = 1<<20; // in bytes
    const int input_file_bytes_thershold = 1358; // in bytes
    const int SWTFC_pool_size_base = 200;
    const int SWTFC_window_size_base = 150;
    std::string input_file_buffer;
    char input_file_chunk[input_file_read_chunk_size];
    std::ifstream input_file(script_input_file, std::ios::binary);
    if (!input_file.is_open()) {
        LOG(ERROR) << "Failed to open input file " << script_input_file;
        return 1;
    }

    const int line_pool_size = 400;
    std::vector <std::string> valid_line_pool;
    valid_line_pool.reserve(line_pool_size);
    std::vector <ULong64_t> timestamp_pool;
    timestamp_pool.reserve(line_pool_size);
    std::vector <UShort_t> fpga_id_pool;
    fpga_id_pool.reserve(line_pool_size);
    std::vector <UShort_t> half_id_pool;
    half_id_pool.reserve(line_pool_size);
    std::vector <bool> line_matched_flags;
    line_matched_flags.reserve(line_pool_size);

    bool flag_continue_reading = true;
    bool input_file_header_found = false;
    bool input_file_header_end_found = false;
    bool input_file_is_last_chunk = false;

    int counter_heartbeat_packet = 0;
    int counter_header_line = 0;
    int counter_complete_20_lines = 0;
    int counter_not_complete_20_lines = 0;
    int counter_not_complete_20_lines_total = 0;
    int counter_invalid_line = 0;
    int counter_valid_line = 0;

    const unsigned char payload_header_byte_0 = 0xAA;
    const unsigned char payload_header_byte_1 = 0x5A;

    UInt_t current_hear_beat_value = 0;

    while ((input_file.read(input_file_chunk, input_file_read_chunk_size) || input_file.gcount() > 0) && flag_continue_reading) {
        input_file_buffer.append(input_file_chunk, input_file.gcount());
        input_file_is_last_chunk = input_file.eof(); // Check if it is the last chunk
        if (input_file_is_last_chunk) {
            LOG(INFO) << "Last chunk found!";
        }
        while ((input_file_buffer.size() > input_file_bytes_thershold || (input_file_is_last_chunk && input_file_buffer.size() >= input_file_bytes_thershold)) && flag_continue_reading) {
            if (input_file_buffer.at(0) == DATA_PACKET_HEADER_START && !(input_file_header_found && input_file_header_end_found)) {
                if (!input_file_header_found) {
                    if (input_file_buffer.at(10) == DATA_PACKET_HEADER_END) {
                        input_file_header_found = true;
                        LOG(INFO) << "Header start found!";
                        size_t pos = input_file_buffer.find('\n');
                        if (pos != std::string::npos) {
                            input_file_buffer.erase(0, pos + 1);
                        }
                        counter_header_line++;
                        continue;
                    }
                }
                if (input_file_header_found && !input_file_header_end_found) {
                    if (input_file_buffer.at(10) == DATA_PACKET_HEADER_END) {
                        input_file_header_end_found = true;
                        LOG(INFO) << "Header end found!";
                        size_t pos = input_file_buffer.find('\n');
                        if (pos != std::string::npos) {
                            input_file_buffer.erase(0, pos + 1);
                        }
                        counter_header_line++;
                        continue;
                    } else {
                        size_t pos = input_file_buffer.find('\n');
                        if (pos == std::string::npos) {
                            break;
                        }
                        counter_header_line++;
                        std::string line = input_file_buffer.substr(0, pos);
                        input_file_buffer.erase(0, pos + 1); // Remove processed line
                        // print out the line in
                        LOG(INFO) << line;
                    }
                }
            } else {
                if (input_file_header_found && input_file_header_end_found) {
                    if (input_file_buffer.at(0) == HEARTBEAT_HEADER_START && input_file_buffer.at(10) == HEARTBEAT_HEADER_END) {
                        input_file_buffer.erase(0, 1358); // Adjust the size to remove the processed frame
                        counter_heartbeat_packet++;
                    } else {
                         // ! -- data packet processing --------------------------------
                        // remove the first 14
                        int _index_base = 0;
                        int _index_end = 1358;
                        // if (input_file_is_last_chunk) {
                        //     if (input_file_buffer.size() <= 1358) {
                        //         _index_end = input_file_buffer.size() - 192;
                        //     }
                        // }
                        // search for the first 0xAA 0x5A header in the first 1358 bytes
                        while (_index_base + 192 <= _index_end) {
                            if ((unsigned char)input_file_buffer.at(_index_base) == payload_header_byte_0 && (unsigned char)input_file_buffer.at(_index_base + 1) == payload_header_byte_1) {
                                break;
                            } else {
                                _index_base++;
                            }
                        }
                        // LOG(INFO) << "First header found at index: " << _index_base;
                        const int frame_size = 192;
                        const int lines_per_event = 4;
                        const int pool_window_size = lines_per_event * 4;
                        while (_index_base + frame_size <= _index_end) {
                            
                            // move the pointer
                            
                            

                            // read one line of 192 bytes
                            std::string line_data(input_file_buffer.begin() + _index_base, input_file_buffer.begin() + _index_base + frame_size);

                            int fpga_id = ((unsigned char) line_data.at(2) & 0xF0) >> 4;
                            
                            int asic_id = (unsigned char) line_data.at(2) & 0x0F;
                            int half_id = (unsigned char) line_data.at(3);
                            // check if the byte# 0 is 0xAA and byte# 1 is 0x5A
                            if ((unsigned char) line_data.at(0) == payload_header_byte_0 && (unsigned char) line_data.at(1) == payload_header_byte_1 && fpga_id <= 0) {
                                counter_valid_line++;
                            } else {
                                if (true) {
                                    counter_invalid_line++;
                                    std::string line_data_hex = "";
                                    std::ostringstream oss;
                                    oss << std::hex << std::setfill('0');
                                    for (int i = 0; i < 20 && i < static_cast<int>(line_data.size()); i++) {
                                        oss << std::setw(2) << static_cast<int>(static_cast<unsigned char>(line_data.at(i))) << " ";
                                    }
                                    line_data_hex = oss.str();
                                    // LOG(WARNING) << "Invalid line found: " << line_data_hex;
                                    _index_base += frame_size;
                                    continue;
                                }
                                // counter_invalid_line++;
                                // // print the first 20 bytes of the invalid line in hex
                                // std::string line_data_hex = "";
                                // std::ostringstream oss;
                                // oss << std::hex << std::setfill('0');
                                // for (int i = 0; i < 20 && i < static_cast<int>(line_data.size()); i++) {
                                //     oss << std::setw(2) << static_cast<int>(static_cast<unsigned char>(line_data.at(i))) << " ";
                                // }
                                // line_data_hex = oss.str();
                                // LOG(WARNING) << "Invalid line found: " << line_data_hex;
                                // continue;
                            }

                            constexpr int timestamp_offset = 16;
                            ULong64_t timestamp = 0;
                            if (line_data.size() < static_cast<size_t>(timestamp_offset + 8)) {
                                LOG(WARNING) << "Line too short for timestamp: size=" << line_data.size();
                                timestamp = 0;
                            } else {
                                for (int i = 0; i < 8; i++) {
                                    timestamp = (timestamp << 8) |
                                                static_cast<ULong64_t>(static_cast<unsigned char>(line_data.at(timestamp_offset + i)));
                                }
                            }
                            // LOG(INFO) << "Frame header: fpga_id=" << fpga_id << " asic_id=" << std::hex << asic_id << " half_id=" << half_id << std::dec << " timestamp=" << timestamp;
                            if (std::find(legal_fpga_id_list.begin(), legal_fpga_id_list.end(), fpga_id) == legal_fpga_id_list.end()) {
                                legal_fpga_id_list.push_back(fpga_id);
                                LOG(INFO) << "New FPGA ID found: " << std::hex << fpga_id << std::dec;
                            }
                            valid_line_pool.push_back(line_data);
                            timestamp_pool.push_back(timestamp);
                            fpga_id_pool.push_back(fpga_id);
                            half_id_pool.push_back(half_id);
                            line_matched_flags.push_back(false);
                            _index_base += frame_size;
                        }
                        input_file_buffer.erase(0, _index_base); // Remove processed data
                        if (valid_line_pool.size() >= static_cast<size_t>(lines_per_event)) {
                            for (int seed_index = 0; seed_index + lines_per_event <= static_cast<int>(valid_line_pool.size()); seed_index++) {
                                if (line_matched_flags.at(seed_index)) {
                                    continue;
                                }
                                auto timestamp_ref = timestamp_pool.at(seed_index);
                                auto fpga_id_ref = fpga_id_pool.at(seed_index);

                                std::vector<int> match_line_index;
                                match_line_index.reserve(lines_per_event);
                                match_line_index.push_back(seed_index);
                                line_matched_flags.at(seed_index) = true;

                                int match_searching_end = seed_index + pool_window_size;
                                if (match_searching_end > static_cast<int>(valid_line_pool.size())) {
                                    match_searching_end = static_cast<int>(valid_line_pool.size());
                                }
                                for (int search_index = seed_index + 1; search_index < match_searching_end; search_index++) {
                                    if (line_matched_flags.at(search_index)) {
                                        continue;
                                    }
                                    if (timestamp_pool.at(search_index) == timestamp_ref && fpga_id_pool.at(search_index) == fpga_id_ref) {
                                        match_line_index.push_back(search_index);
                                        line_matched_flags.at(search_index) = true;
                                        if (static_cast<int>(match_line_index.size()) == lines_per_event) {
                                            break;
                                        }
                                    }
                                }
                                int matched_frame_count = static_cast<int>(match_line_index.size());
                                if (matched_frame_count == lines_per_event) {
                                    counter_complete_20_lines += matched_frame_count;
                                    if (counter_complete_20_lines / lines_per_event >= script_n_events && script_n_events > 0) {
                                        flag_continue_reading = false;
                                        break;
                                    }
                                    std::fill(branch_daqh_list, branch_daqh_list + 4, 0u);
                                    std::fill(branch_crc32_list, branch_crc32_list + 4, 0u);
                                    std::fill(branch_tc_list, branch_tc_list + FPGA_CHANNEL_NUMBER, false);
                                    std::fill(branch_tp_list, branch_tp_list + FPGA_CHANNEL_NUMBER, false);
                                    std::fill(branch_val0_list, branch_val0_list + FPGA_CHANNEL_NUMBER, 0u);
                                    std::fill(branch_val1_list, branch_val1_list + FPGA_CHANNEL_NUMBER, 0u);
                                    std::fill(branch_val2_list, branch_val2_list + FPGA_CHANNEL_NUMBER, 0u);
                                    for (auto &index : match_line_index) {
                                        const auto &matched_line = valid_line_pool.at(index);
                                        int fpga_id = ((unsigned char) matched_line.at(2) & 0xF0) >> 4;
                                        int asic_id = (unsigned char) matched_line.at(2) & 0x0F;
                                        int half_id = (unsigned char) matched_line.at(3);
                                        int half_index = half_id - 0x24;
                                        if (half_index < 0 || half_index >= 4) {
                                            continue;
                                        }
                                        int half_channel_offset = half_index * 38;

                                        if (asic_id == 0x00 && half_id == 0x24) {
                                            *branch_fpga_id = fpga_id;
                                            *branch_timestamp = timestamp_ref;
                                            *branch_last_heartbeat = current_hear_beat_value;
                                        }

                                        // handle the data line with 160 bytes of payload
                                        constexpr int payload_offset = 32;
                                        constexpr int payload_size = 160;
                                        constexpr int payload_words = payload_size / 4;
                                        if (matched_line.size() < static_cast<size_t>(payload_offset + payload_size)) {
                                            continue;
                                        }

                                        auto read_word = [&](int word_index) {
                                            int base = payload_offset + word_index * 4;
                                            return (static_cast<ULong64_t>(static_cast<unsigned char>(matched_line.at(base))) << 24) |
                                                   (static_cast<ULong64_t>(static_cast<unsigned char>(matched_line.at(base + 1))) << 16) |
                                                   (static_cast<ULong64_t>(static_cast<unsigned char>(matched_line.at(base + 2))) << 8) |
                                                   static_cast<ULong64_t>(static_cast<unsigned char>(matched_line.at(base + 3)));
                                        };

                                        UInt_t daqh = static_cast<UInt_t>(read_word(0));
                                        branch_daqh_list[half_index] = daqh;

                                        for (int word_index = 1; word_index <= 38; word_index++) {
                                            ULong64_t word = read_word(word_index);
                                            Bool_t tc = (word >> 31) & 0x1;
                                            Bool_t tp = (word >> 30) & 0x1;
                                            UInt_t val0 = (word >> 20) & 0x3ff;
                                            UInt_t val1 = (word >> 10) & 0x3ff;
                                            UInt_t val2 = word & 0x3ff;
                                            int channel_index = word_index - 1;
                                            branch_tc_list[channel_index + half_channel_offset] = tc;
                                            branch_tp_list[channel_index + half_channel_offset] = tp;
                                            branch_val0_list[channel_index + half_channel_offset] = val0;
                                            branch_val1_list[channel_index + half_channel_offset] = val1;
                                            branch_val2_list[channel_index + half_channel_offset] = val2;
                                        }

                                        UInt_t crc32 = static_cast<UInt_t>(read_word(payload_words - 1));
                                        branch_crc32_list[half_index] = crc32;
                                    }
                                    output_tree->Fill();
                                } else {
                                    int matched_frame_count = static_cast<int>(match_line_index.size());
                                    counter_not_complete_20_lines += matched_frame_count;
                                    counter_not_complete_20_lines_total += 1;
                                }
                            }
                            if (static_cast<int>(valid_line_pool.size()) > pool_window_size) {
                                int erase_count = static_cast<int>(valid_line_pool.size()) - pool_window_size;
                                valid_line_pool.erase(valid_line_pool.begin(), valid_line_pool.begin() + erase_count);
                                timestamp_pool.erase(timestamp_pool.begin(), timestamp_pool.begin() + erase_count);
                                fpga_id_pool.erase(fpga_id_pool.begin(), fpga_id_pool.begin() + erase_count);
                                half_id_pool.erase(half_id_pool.begin(), half_id_pool.begin() + erase_count);
                                line_matched_flags.erase(line_matched_flags.begin(), line_matched_flags.begin() + erase_count);
                            }
                        } // end of if (valid_line_pool.size() >= lines_per_event)
                    }
                }
            }
        }
    }

    LOG(INFO) << "Total header lines: " << counter_header_line;
    LOG(INFO) << "Total heartbeat packets: " << counter_heartbeat_packet;
    LOG(INFO) << "Total valid lines: " << counter_valid_line;
    LOG(INFO) << "Total invalid lines: " << counter_invalid_line;
    if (legal_fpga_id_list.size() == 0) {
        LOG(ERROR) << "No legal FPGA ID found!";
    } else {
        LOG(INFO) << "Total complete 20 lines: " << counter_complete_20_lines << " (" << counter_complete_20_lines / legal_fpga_id_list.size() / 20 << " samples per FPGA)";
        LOG(INFO) << "Total not complete 20 lines: " << counter_not_complete_20_lines << " (" << counter_not_complete_20_lines_total / legal_fpga_id_list.size() << " samples per FPGA)";
    }
    if ((counter_complete_20_lines + counter_not_complete_20_lines) == 0) {
        LOG(ERROR) << "No valid line found!";
    } else {
        LOG(INFO) << "Complete 20 lines percentage: " << (float)counter_complete_20_lines / (counter_complete_20_lines + counter_not_complete_20_lines) * 100 << "%";
        LOG(INFO) << "Total valid line percentage: " << (float)counter_valid_line / (counter_valid_line + counter_invalid_line / 40) * 100 << "%";
    }
    std::string _legal_fpga_id_list_str = "";
    for (auto &fpga_id : legal_fpga_id_list) {
        _legal_fpga_id_list_str += std::to_string(fpga_id) + " ";
    }
    LOG(INFO) << "Legal FPGA ID list: " << _legal_fpga_id_list_str;

    input_file.close();

    // * --- Write output file ----------------------------------------------------------
    // * --------------------------------------------------------------------------------
    // -- Write the tree
    output_tree->Write();
    // -- Write the meta data
    TNamed("Rootifier_script_name", script_name.c_str()).Write();
    TNamed("Rootifier_script_version", script_version.c_str()).Write();
    TNamed("Rootifier_script_input_file", script_input_file.c_str()).Write();
    TNamed("Rootifier_script_output_file", script_output_file.c_str()).Write();
    TNamed("Rootifier_script_n_events", std::to_string(script_n_events).c_str()).Write();
    TNamed("Rootifier_script_output_folder", script_output_folder.c_str()).Write();
    // -- Write the runnning info
    TNamed("Rootifier_counter_header_line", std::to_string(counter_header_line).c_str()).Write();
    TNamed("Rootifier_counter_heartbeat_packet", std::to_string(counter_heartbeat_packet).c_str()).Write();
    TNamed("Rootifier_counter_complete_20_lines", std::to_string(counter_complete_20_lines).c_str()).Write();
    TNamed("Rootifier_counter_not_complete_20_lines", std::to_string(counter_not_complete_20_lines).c_str()).Write();
    TNamed("Rootifier_counter_not_complete_20_lines_total", std::to_string(counter_not_complete_20_lines_total).c_str()).Write();
    TNamed("Rootifier_counter_invalid_line", std::to_string(counter_invalid_line).c_str()).Write();
    TNamed("Rootifier_counter_valid_line", std::to_string(counter_valid_line).c_str()).Write();
    // -- Write the legal fpga id list
    TNamed("Rootifier_legal_fpga_id_list", _legal_fpga_id_list_str.c_str()).Write();
    
    output_root->Close();
    return 0;
}