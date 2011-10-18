#include "held_karp.h"
#include "utils.h"

#include "mongoose.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <cstdio>
#include <cstring>

const int MIN_N = 10;
const int MAX_N = 300;

std::string processRequest(int n) {
    std::vector<std::pair<double, double> > points = tsp::utils::generateEuclidean(n);
    std::vector<std::vector<double> > metric = tsp::utils::getEuclideanMetric(points);
    std::vector<std::vector<char> > partialTour(n, std::vector<char>(n, -1));
    tsp::core::HeldKarpLowerBound solver(metric, partialTour);
    std::vector<std::vector<double> > fractionalTour(n, std::vector<double>(n));
    std::cerr << "computing Held-Karp... " << std::endl;
    solver.getBound(fractionalTour);
    std::cerr << "done" << std::endl;
    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss.precision(5);
    oss << "{";
    oss << "\"points\": [";
    for (int i = 0; i < n; ++i) {
        if (i != 0) {
            oss << ",";
        }
        oss << "[" << points[i].first << "," << points[i].second << "]";
    }
    oss << "],";
    oss << "\"fixedEdges\": [";
    bool first = true;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (fractionalTour[i][j] > 1.0 - tsp::utils::EPSILON) {
                if (first) {
                    first = false;
                }
                else {
                    oss << ",";
                }
                oss << "[" << i << "," << j << "]";
            }
        }
    }
    oss << "],";
    oss << "\"vagueEdges\": [";
    first = true;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double value = fractionalTour[i][j];
            if (value > tsp::utils::EPSILON && value < 1.0 - tsp::utils::EPSILON) {
                if (first) {
                    first = false;
                }
                else {
                    oss << ",";
                }
                oss << "[" << i << "," << j << "," << value << "]";
            }
        }
    }
    oss << "]";
    oss << "}";
    return oss.str();
}

void sendFile(struct mg_connection *conn, const std::string &fileName) {
    FILE *input = fopen(fileName.c_str(), "rb");
    fseek(input, 0, SEEK_END);
    size_t fileSize = ftell(input);
    rewind(input);
    char *buffer = new char[fileSize];
    fread(buffer, 1, fileSize, input);
    fclose(input);
    mg_write(conn, buffer, fileSize);
    delete[] buffer;
}

void *callback(enum mg_event event,
               struct mg_connection *conn,
               const struct mg_request_info *request_info) {
    if (event == MG_NEW_REQUEST) {
        std::string uri(request_info->uri);
        std::cerr << "uri: " << uri << std::endl;
        if (uri == "/") {
            mg_printf(conn, "HTTP/1.1 200 OK\nContent-Type: text/html\n\n");
            sendFile(conn, "frontend.html");
            return (void*)1;
        }
        else if (uri == "/engine") {
            if (!request_info->query_string) {
                return NULL;
            }
            mg_printf(conn, "HTTP/1.1 200 OK\nContent-Type: text/plain\n\n");
            char buf[100];
            mg_get_var(request_info->query_string, strlen(request_info->query_string), "n", buf, sizeof(buf));
            if (buf[0] == '\0') {
                return NULL;
            }
            std::stringstream ss(buf);
            int n;
            if (!(ss >> n)) {
                return NULL;
            }
            if (n < MIN_N && n > MAX_N) {
                return NULL;
            }
            std::cerr << "good n: " << n << std::endl;
            std::string response = processRequest(n);
            mg_write(conn, response.c_str(), response.size());
            return (void*)1;
        }
        return NULL;
    }
    else {
        return NULL;
    }
}

int main() {
  struct mg_context *ctx;
  const char *options[] = {"listening_ports", "8080", NULL};

  ctx = mg_start(&callback, NULL, options);
  getchar();
  mg_stop(ctx);

  return 0;
}
