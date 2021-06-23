/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   IP, Inti Pelupessy,  NL eScience Center, i.pelupessy@esciencecenter.nl  *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#include <MultiProcess/BidirMMapPipe.h>

#include <sys/time.h>
#include <iomanip>

#include "gtest/gtest.h"

using namespace RooFit;

int simplechild(BidirMMapPipe& pipe)
{
    // child does an echo loop
    while (pipe.good() && !pipe.eof()) {
        // read a string
        std::string str;
        pipe >> str;
        if (!pipe) return -1;
        if (pipe.eof()) break;
        if (!str.empty()) {
            std::cout << "[CHILD (PID " << getpid() << ")] :  read: " << str << std::endl;
            str = str + "... early in the morning?";
        }
        pipe << str << BidirMMapPipe::flush;
        // did our parent tell us to shut down?
        if (str.empty()) break;
        if (!pipe) return -1;
        if (pipe.eof()) break;
        std::cout << "[CHILD (PID " << getpid() << ")] : wrote: " << str << std::endl;
    }
    std::cout << "[CHILD (PID " << getpid() << ")] : shutting down "  << pipe.isParent() << std::endl;
    std::cout << "[CHILD (PID " << getpid() << ")] : close " << pipe.close() << std::endl;
    return 0;
}


#include <sstream>
int randomchild(BidirMMapPipe& pipe, bool grand)
{
    // child sends out something at random intervals
    ::srand48(::getpid());
    {
        // wait for parent's go ahead signal
        std::string s;
        pipe >> s;
    }
    std::string prestring, PRESTRING;
    if (grand) {
        prestring = "grand";
        PRESTRING = "GRAND";
    }
    // no shutdown sequence needed on this side - we're producing the data,
    // and the parent can just read until we're done (when it'll get EOF)
    for (int i = 0; i < 5; ++i) {
        // sleep a random time between 0 and .9 seconds
        ::usleep(int(1e6 * ::drand48()));
        std::ostringstream buf;
        buf << prestring << "child pid " << ::getpid() << " sends message " << i;
        std::string str = buf.str();
        std::cout << "[" << PRESTRING << "CHILD (PID " << getpid() << ")] : " << str << std::endl;
        pipe << str << BidirMMapPipe::flush;
        if (!pipe) return -1;
        if (pipe.eof()) break;
    }
    // tell parent we're shutting down
    pipe << "" << BidirMMapPipe::flush;
    // wait for parent to acknowledge
    std::string s;
    pipe >> s;
    pipe.close();
    return 0;
}

int benchchildrtt(BidirMMapPipe& pipe)
{
    // child does the equivalent of listening for pings and sending the
    // packet back
    char* str = 0;
    while (pipe && !pipe.eof()) {
        pipe >> str;
        if (!pipe) {
            std::free(str);
            pipe.close();
            return -1;
        }
        if (pipe.eof()) break;
        pipe << str << BidirMMapPipe::flush;
        // if we have just completed the shutdown handshake, we break here
        if (!std::strlen(str)) break;
    }
    std::free(str);
    pipe.close();
    return 0;
}

int benchchildsink(BidirMMapPipe& pipe)
{
    // child behaves like a sink
    char* str = 0;
    while (pipe && !pipe.eof()) {
        pipe >> str;
        if (!std::strlen(str)) break;
    }
    pipe << "" << BidirMMapPipe::flush;
    std::free(str);
    pipe.close();
    return 0;
}

int benchchildsource(BidirMMapPipe& pipe)
{
    // child behaves like a source
    char* str = 0;
    for (unsigned i = 0; i <= 24; ++i) {
        str = reinterpret_cast<char*>(std::realloc(str, (1 << i) + 1));
        std::memset(str, '4', 1 << i);
        str[1 << i] = 0;
        for (unsigned j = 0; j < 1 << 7; ++j) {
            pipe << str;
            if (!pipe || pipe.eof()) {
                std::free(str);
                pipe.close();
                return -1;
            }
        }
        // tell parent we're done with this block size
        pipe << "" << BidirMMapPipe::flush;
    }
    // tell parent to shut down
    pipe << "" << BidirMMapPipe::flush;
    std::free(str);
    pipe.close();
    return 0;
}

template <typename... Args>
BidirMMapPipe* spawnChild(int (*childexec)(BidirMMapPipe&, Args...), Args... args)
{
    // create a pipe with the given child at the remote end
    BidirMMapPipe *p = new BidirMMapPipe();
    if (p->isChild()) {
        int retVal = childexec(*p, args...);
        delete p;
        std::_Exit(retVal);
    }
    return p;
}

// simple echo loop test
TEST(BidirMMapPipe, simple){
    std::cout << "[PARENT]: simple challenge-response test, "
        "one child:" << std::endl;
    BidirMMapPipe* pipe = spawnChild(simplechild);
    for (int i = 0; i < 5; ++i) {
        std::string str("What shall we do with a drunken sailor...");
        *pipe << str << BidirMMapPipe::flush;
        ASSERT_TRUE(*pipe);
        std::cout << "[PARENT (PID " << getpid() << ")]: wrote: " << str << std::endl;
        *pipe >> str;
        ASSERT_TRUE(*pipe);
        std::cout << "[PARENT (PID " << getpid() << ")]:  read: " << str << std::endl;
    }
    // send shutdown string
    *pipe << "" << BidirMMapPipe::flush;
    // wait for shutdown handshake
    std::string s;
    *pipe >> s;
    int retVal = pipe->close();
    std::cout << "[PARENT (PID " << getpid() << ")]: exit status of child: " << retVal <<
        std::endl;
    delete pipe;
}

    // simple poll test - children send 5 results in random intervals
TEST(BidirMMapPipe, poll)
    {
        unsigned nch = 20;
        std::cout << std::endl << "[PARENT (PID " << getpid() << ")]: polling test, " << nch <<
            " children:" << std::endl;
        typedef BidirMMapPipe::PollEntry PollEntry;
        // poll data structure
        BidirMMapPipe::PollVector pipes;
        pipes.reserve(nch);
        // spawn children
        for (unsigned i = 0; i < nch; ++i) {
            std::cout << "[PARENT (PID " << getpid() << ")]: spawning child " << i << std::endl;
            pipes.push_back(PollEntry(spawnChild(randomchild, false),
                        BidirMMapPipe::Readable));
        }
        // wake children up
        std::cout << "[PARENT (PID " << getpid() << ")]: waking up children" << std::endl;
        for (unsigned i = 0; i < nch; ++i)
            *pipes[i].pipe << "" << BidirMMapPipe::flush;
        std::cout << "[PARENT (PID " << getpid() << ")]: waiting for events on children's pipes" << std::endl;
        // while at least some children alive
        while (!pipes.empty()) {
            // poll, wait until status change (infinite timeout)
            int npipes = BidirMMapPipe::poll(pipes, -1);
            // scan for pipes with changed status
            for (std::vector<PollEntry>::iterator it = pipes.begin();
                    npipes && pipes.end() != it; ) {
                if (!it->revents) {
                    // unchanged, next one
                    ++it;
                    continue;
                }
                --npipes; // maybe we can stop early...
                // read from pipes which are readable
                if (it->revents & BidirMMapPipe::Readable) {
                    std::string s;
                    *(it->pipe) >> s;
                    if (!s.empty()) {
                        std::cout << "[PARENT (PID " << getpid() << ")]: Read from pipe " << it->pipe <<
                            ": " << s << std::endl;
                        ++it;
                        continue;
                    } else {
                        // child is shutting down...
                        *(it->pipe) << "" << BidirMMapPipe::flush;
                        goto childcloses;
                    }
                }
                // retire pipes with error or end-of-file condition
                if (it->revents & (BidirMMapPipe::Error |
                            BidirMMapPipe::EndOfFile |
                            BidirMMapPipe::Invalid)) {
                    std::cerr << "[DEBUG]: Event on pipe " << it->pipe <<
                        " revents" <<
                        ((it->revents & BidirMMapPipe::Readable) ? " Readable" : "") <<
                        ((it->revents & BidirMMapPipe::Writable) ? " Writable" : "") <<
                        ((it->revents & BidirMMapPipe::ReadError) ? " ReadError" : "") <<
                        ((it->revents & BidirMMapPipe::WriteError) ? " WriteError" : "") <<
                        ((it->revents & BidirMMapPipe::ReadEndOfFile) ? " ReadEndOfFile" : "") <<
                        ((it->revents & BidirMMapPipe::WriteEndOfFile) ? " WriteEndOfFile" : "") <<
                        ((it->revents & BidirMMapPipe::ReadInvalid) ? " ReadInvalid" : "") <<
                        ((it->revents & BidirMMapPipe::WriteInvalid) ? " WriteInvalid" : "") <<
                        std::endl;
childcloses:
                    int retVal = it->pipe->close();
                    std::cout << "[PARENT (PID " << getpid() << ")]: child exit status: " <<
                        retVal << ", number of children still alive: " <<
                        (pipes.size() - 1) << std::endl;
                    delete it->pipe;
                    it = pipes.erase(it);
                    continue;
                }
            }
        }
    }
    //~ // little benchmark - round trip time
    //~ {
        //~ std::cout << std::endl << "[PARENT]: benchmark: round-trip times vs block size" << std::endl;
        //~ for (unsigned i = 0; i <= 24; ++i) {
            //~ char *s = new char[1 + (1 << i)];
            //~ std::memset(s, 'A', 1 << i);
            //~ s[1 << i] = 0;
            //~ const unsigned n = 1 << 7;
            //~ double avg = 0., min = 1e42, max = -1e42;
            //~ BidirMMapPipe *pipe = spawnChild(benchchildrtt);
            //~ for (unsigned j = n; j--; ) {
                //~ struct timeval t1;
                //~ ::gettimeofday(&t1, 0);
                //~ *pipe << s << BidirMMapPipe::flush;
                //~ if (!*pipe || pipe->eof()) break;
                //~ *pipe >> s;
                //~ if (!*pipe || pipe->eof()) break;
                //~ struct timeval t2;
                //~ ::gettimeofday(&t2, 0);
                //~ t2.tv_sec -= t1.tv_sec;
                //~ t2.tv_usec -= t1.tv_usec;
                //~ double dt = 1e-6 * double(t2.tv_usec) + double(t2.tv_sec);
                //~ if (dt < min) min = dt;
                //~ if (dt > max) max = dt;
                //~ avg += dt;
            //~ }
            //~ // send a shutdown string
            //~ *pipe << "" << BidirMMapPipe::flush;
            //~ // get child's shutdown ok
            //~ *pipe >> s;
            //~ avg /= double(n);
            //~ avg *= 1e6; min *= 1e6; max *= 1e6;
            //~ int retVal = pipe->close();
            //~ if (retVal) {
                //~ std::cout << "[PARENT]: child exited with code " << retVal << std::endl;
                //~ return retVal;
            //~ }
            //~ delete pipe;
            //~ // there is a factor 2 in the formula for the transfer rate below,
            //~ // because we transfer data of twice the size of the block - once
            //~ // to the child, and once for the return trip
            //~ std::cout << "block size " << std::setw(9) << (1 << i) <<
                //~ " avg " << std::setw(7) << avg << " us min " <<
                //~ std::setw(7) << min << " us max " << std::setw(7) << max <<
                //~ "us speed " << std::setw(9) <<
                //~ 2. * (double(1 << i) / double(1 << 20) / (1e-6 * avg)) <<
                //~ " MB/s" << std::endl;
            //~ delete[] s;
        //~ }
        //~ std::cout << "[PARENT]: all children had exit code 0" << std::endl;
    //~ }
    //~ // little benchmark - child as sink
    //~ {
        //~ std::cout << std::endl << "[PARENT]: benchmark: raw transfer rate with child as sink" << std::endl;
        //~ for (unsigned i = 0; i <= 24; ++i) {
            //~ char *s = new char[1 + (1 << i)];
            //~ std::memset(s, 'A', 1 << i);
            //~ s[1 << i] = 0;
            //~ const unsigned n = 1 << 7;
            //~ double avg = 0., min = 1e42, max = -1e42;
            //~ BidirMMapPipe *pipe = spawnChild(benchchildsink);
            //~ for (unsigned j = n; j--; ) {
                //~ struct timeval t1;
                //~ ::gettimeofday(&t1, 0);
                //~ // streaming mode - we do not flush here
                //~ *pipe << s;
                //~ if (!*pipe || pipe->eof()) break;
                //~ struct timeval t2;
                //~ ::gettimeofday(&t2, 0);
                //~ t2.tv_sec -= t1.tv_sec;
                //~ t2.tv_usec -= t1.tv_usec;
                //~ double dt = 1e-6 * double(t2.tv_usec) + double(t2.tv_sec);
                //~ if (dt < min) min = dt;
                //~ if (dt > max) max = dt;
                //~ avg += dt;
            //~ }
            //~ // send a shutdown string
            //~ *pipe << "" << BidirMMapPipe::flush;
            //~ // get child's shutdown ok
            //~ *pipe >> s;
            //~ avg /= double(n);
            //~ avg *= 1e6; min *= 1e6; max *= 1e6;
            //~ int retVal = pipe->close();
            //~ if (retVal) {
                //~ std::cout << "[PARENT]: child exited with code " << retVal << std::endl;
                //~ return retVal;
            //~ }
            //~ delete pipe;
            //~ std::cout << "block size " << std::setw(9) << (1 << i) <<
                //~ " avg " << std::setw(7) << avg << " us min " <<
                //~ std::setw(7) << min << " us max " << std::setw(7) << max <<
                //~ "us speed " << std::setw(9) <<
                //~ (double(1 << i) / double(1 << 20) / (1e-6 * avg)) <<
                //~ " MB/s" << std::endl;
            //~ delete[] s;
        //~ }
        //~ std::cout << "[PARENT]: all children had exit code 0" << std::endl;
    //~ }
    //~ // little benchmark - child as source
    //~ {
        //~ std::cout << std::endl << "[PARENT]: benchmark: raw transfer rate with child as source" << std::endl;
        //~ char *s = 0;
        //~ double avg = 0., min = 1e42, max = -1e42;
        //~ unsigned n = 0, bsz = 0;
        //~ BidirMMapPipe *pipe = spawnChild(benchchildsource);
        //~ while (*pipe && !pipe->eof()) {
            //~ struct timeval t1;
            //~ ::gettimeofday(&t1, 0);
            //~ // streaming mode - we do not flush here
            //~ *pipe >> s;
            //~ if (!*pipe || pipe->eof()) break;
            //~ struct timeval t2;
            //~ ::gettimeofday(&t2, 0);
            //~ t2.tv_sec -= t1.tv_sec;
            //~ t2.tv_usec -= t1.tv_usec;
            //~ double dt = 1e-6 * double(t2.tv_usec) + double(t2.tv_sec);
            //~ if (std::strlen(s)) {
                //~ ++n;
                //~ if (dt < min) min = dt;
                //~ if (dt > max) max = dt;
                //~ avg += dt;
                //~ bsz = std::strlen(s);
            //~ } else {
                //~ if (!n) break;
                //~ // next block size
                //~ avg /= double(n);
                //~ avg *= 1e6; min *= 1e6; max *= 1e6;
//~ 
                //~ std::cout << "block size " << std::setw(9) << bsz <<
                    //~ " avg " << std::setw(7) << avg << " us min " <<
                    //~ std::setw(7) << min << " us max " << std::setw(7) <<
                    //~ max << "us speed " << std::setw(9) <<
                    //~ (double(bsz) / double(1 << 20) / (1e-6 * avg)) <<
                    //~ " MB/s" << std::endl;
                //~ n = 0;
                //~ avg = 0.;
                //~ min = 1e42;
                //~ max = -1e42;
            //~ }
        //~ }
        //~ int retVal = pipe->close();
            //~ std::cout << "[PARENT]: child exited with code " << retVal << std::endl;
        //~ if (retVal) return retVal;
        //~ delete pipe;
        //~ std::free(s);
    //~ }
    //~ return 0;

// additional tests


int simplerelay(BidirMMapPipe& in,BidirMMapPipe& out)
{
    // child does an echo loop
    while (in.good() && !in.eof()) {
        // read a string
        std::string str;
        in >> str;
        if (!in) return -1;
        if (in.eof()) break;
        //~ if (!str.empty()) {
            std::cout << "[RELAY (PID " << getpid() << ")] :  from in: " << str << std::endl;
            out << str << BidirMMapPipe::flush;
            std::cout << "[RELAY (PID " << getpid() << ")] : to out: " << str << std::endl;
        //~ }
        // did our parent tell us to shut down?
        if (str.empty()) break;
        if (!out) return -1;
        if (out.eof()) break;
        out >> str;
        if (!in) return -1;
        if (in.eof()) break;
        //~ if (!str.empty()) {
            std::cout << "[RELAY (PID " << getpid() << ")] :  from out: " << str << std::endl;
            in << str << BidirMMapPipe::flush;
            std::cout << "[RELAY (PID " << getpid() << ")] :  to in: " << str << std::endl;
        //~ }
    }
    std::cout << "[RELAY (PID " << getpid() << ")] : shutting down "  << std::endl;
    return 0;
}


#include <sys/time.h>
#include <iomanip>

// simple echo loop test
template<class T>
int simple_echo(T& out, BidirMMapPipe* pipe, std::string extra_string = "")
{
    {
        out << "[PARENT (PID " << getpid() << ")]: simple challenge-response test, "
            "one child (extra_string: " << extra_string << "):" << std::endl;
        for (int i = 0; i < 2; ++i) {
            std::string str("What shall we do with a drunken sailor...");
            str += extra_string;
            *pipe << str + "  "+ std::to_string(i) << BidirMMapPipe::flush;
            if (!*pipe) return -1;
            out << "[PARENT (PID " << getpid() << ")]: wrote: " << str << std::endl;
            *pipe >> str;
            if (!*pipe) return -1;
            out << "[PARENT (PID " << getpid() << ")]:  read: " << str << std::endl;
//            out.flush();
        }
        // send shutdown string
        *pipe << "" << BidirMMapPipe::flush;
        // wait for shutdown handshake
        std::string s;
        *pipe >> s;
        int retVal = pipe->close();
        out << "[PARENT (PID " << getpid() << ")]: status of child: " << std::to_string(retVal) <<
                                                                                                std::endl;
//        out.flush();
        return retVal;
    }  
}

void simple_echo_direct()
{
    BidirMMapPipe* pipe = spawnChild(simplechild);
    simple_echo(std::cout, pipe);
    delete pipe;
}

void simple_echo_relay()
{
    BidirMMapPipe* pipe1 = spawnChild(simplechild);
    BidirMMapPipe* pipe2 = new BidirMMapPipe(true, false, false);

    if(pipe2->isChild()) {
      std :: cout << "child of p2 considers itself pipe parent? " << pipe1->isParent() << std:: endl;
      simplerelay(*pipe2, *pipe1); // forward output over pipe2
      pipe2->close();
      pipe1->close();
    }
    if(pipe2->isParent()) {
      simple_echo(std::cout, pipe2);
    }
    std::cout << "ending" << std::endl;
}

TEST(BidirMMapPipe, direct)
{
    simple_echo_direct();
}


TEST(BidirMMapPipe, relay)
{
    simple_echo_relay();
}


TEST(BidirMMapPipe, bothDirectAndRelay)
{
    simple_echo_direct();
    simple_echo_relay();
    simple_echo_direct();
}


TEST(BidirMMapPipe, grandChild) {
    BidirMMapPipe* pipe1 = new BidirMMapPipe();
    if (pipe1->isChild()) {
        BidirMMapPipe* pipe2 = spawnChild(simplechild);

        if(pipe2->isParent()) {
            simplerelay(*pipe1, *pipe2); // forward output over pipe2
            pipe1->close();
        }
    } else {
        // send an echo over pipe1 from master, which simplerelay should send on over pipe2
        simple_echo(std::cout, pipe1);
    }

    std::cout << "ending" << std::endl;
}

TEST(BidirMMapPipe, greatGrandChild) {
    BidirMMapPipe* pipe1 = new BidirMMapPipe();
    if (pipe1->isChild()) {
        BidirMMapPipe* pipe2 = new BidirMMapPipe();

        if(pipe2->isParent()) {
            simplerelay(*pipe1, *pipe2); // forward output over pipe2
        } else {
            BidirMMapPipe *pipe3 = spawnChild(simplechild);

            if (pipe3->isParent()) {
                simplerelay(*pipe2, *pipe3); // forward output over pipe2
                pipe2->close();
                pipe1->close();
            }
        }

    } else {
        // send an echo over pipe1 from master, which simplerelay should send on over pipe2 and pipe3
        simple_echo(std::cout, pipe1);
    }

    std::cout << "ending" << std::endl;
}


TEST(BidirMMapPipe, multipleChildrenAndGrandChildren) {
    BidirMMapPipe* child1 = new BidirMMapPipe();
    if (child1->isChild()) {
        BidirMMapPipe* grandchild1_1 = spawnChild(simplechild);

        if(grandchild1_1->isParent()) {
            simplerelay(*child1, *grandchild1_1); // forward output over grandchild1_1
            delete grandchild1_1;
            std::_Exit(0); // exit child1 process
        }
    } else {
        BidirMMapPipe* child2 = new BidirMMapPipe();

        if (child2->isChild()) {
            BidirMMapPipe* grandchild2_1 = spawnChild(simplechild);

            if(grandchild2_1->isParent()) {
                simplerelay(*child2, *grandchild2_1); // forward output over grandchild2_1
//                child2->close();
                delete grandchild2_1;
//                delete child2;
                std::_Exit(0); // exit child2 process
            }
        } else {
            // send an echo over child2 from master, which simplerelay should send on over grandchild2_1
            simple_echo(std::cout, child2, " second child");
        }

        delete child2;

        // send an echo over child1 from master, which simplerelay should send on over grandchild1_1
        simple_echo(std::cout, child1, " first child");
    }

    delete child1;

    std::cout << "ending" << std::endl;
}


void poll_relay(BidirMMapPipe& parent, BidirMMapPipe::PollVector & grandchildren) {
    {
        // wait for parent's go ahead signal
        std::string s;
        parent >> s;
        // and relay it to grandchildren
        std::cout << "[CHILD (PID " << getpid() << ")]: waking up my grandchildren" << std::endl;
        for (unsigned i = 0; i < grandchildren.size(); ++i) {
            *grandchildren[i].pipe << "" << BidirMMapPipe::flush;
        }
    }

    std::cout << "[CHILD (PID " << getpid() << ")]: waiting for events on grandchildren's pipes" << std::endl;
    // while at least some children alive
    while (!grandchildren.empty()) {
        // poll, wait until status change (infinite timeout)
        int npipes = BidirMMapPipe::poll(grandchildren, -1);
        // scan for pipes with changed status
        for (auto it = grandchildren.begin(); npipes && grandchildren.end() != it; ) {
            if (!it->revents) {
                // unchanged, next one
                ++it;
                continue;
            }
            --npipes; // maybe we can stop early...
            // read from pipes which are readable
            if (it->revents & BidirMMapPipe::Readable) {
                std::string s;
                pid_t grandchild_pid = it->pipe->pidOtherEnd();
                *(it->pipe) >> s;
                if (!s.empty()) {
                    std::cout << "[CHILD (PID " << getpid() << ")]: Read from pipe " << it->pipe <<
                              ": " << s << ", passing on to parent" << std::endl;
                    std::stringstream ss;
                    ss << "from child " << getpid() << s;
                    parent << ss.str() << grandchild_pid << BidirMMapPipe::flush;
                    ++it;
                    continue;
                } else {
                    // grandchild is shutting down...
                    *(it->pipe) << "" << BidirMMapPipe::flush;
                    goto grandchildcloses;
                }
            }
            // retire pipes with error or end-of-file condition
            if (it->revents & (BidirMMapPipe::Error |
                               BidirMMapPipe::EndOfFile |
                               BidirMMapPipe::Invalid)) {
                std::cerr << "[DEBUG GRANDCHILD]: Event on pipe " << it->pipe <<
                          " revents" <<
                          ((it->revents & BidirMMapPipe::Readable) ? " Readable" : "") <<
                          ((it->revents & BidirMMapPipe::Writable) ? " Writable" : "") <<
                          ((it->revents & BidirMMapPipe::ReadError) ? " ReadError" : "") <<
                          ((it->revents & BidirMMapPipe::WriteError) ? " WriteError" : "") <<
                          ((it->revents & BidirMMapPipe::ReadEndOfFile) ? " ReadEndOfFile" : "") <<
                          ((it->revents & BidirMMapPipe::WriteEndOfFile) ? " WriteEndOfFile" : "") <<
                          ((it->revents & BidirMMapPipe::ReadInvalid) ? " ReadInvalid" : "") <<
                          ((it->revents & BidirMMapPipe::WriteInvalid) ? " WriteInvalid" : "") <<
                          std::endl;
              grandchildcloses:
                int retVal = it->pipe->close();
                std::cout << "[CHILD (PID " << getpid() << ")]: grandchild exit status: " <<
                          retVal << ", number of grandchildren still alive: " <<
                          (grandchildren.size() - 1) << std::endl;
                delete it->pipe;
                it = grandchildren.erase(it);
                continue;
            }
        }
    }

}


// hierarchical poll test - grandchildren send 5 results in random intervals
TEST(BidirMMapPipe, pollHierarchy) {
    unsigned nch = 5;  // spawn 5 children and 5x5 grandchildren, 31 processes in total
    std::cout << std::endl << "[PARENT (PID " << getpid() << ")]: polling test, " << nch <<
              " children, " << nch*nch << " grandchildren:" << std::endl;

    // poll data structures
    BidirMMapPipe::PollVector children, grandchildren;
    children.reserve(nch);
    grandchildren.reserve(nch);
    std::map<pid_t, unsigned> N_received_signals;

    // spawn children
    for (unsigned i = 0; i < nch; ++i) {
        std::cout << "[PARENT (PID " << getpid() << ")]: spawning child " << i << std::endl;
        children.emplace_back(new BidirMMapPipe(), BidirMMapPipe::Readable);

        BidirMMapPipe& child = *children.back().pipe;
        // THE BELOW BLOCK TAKES PLACE ON THE CHILDREN
        if (child.isChild()) {
            for (unsigned j = 0; j < nch; ++j) {
                std::cout << "[CHILD " << i << " (PID " << getpid() << ")]: spawning grandchild " << j << std::endl;
                grandchildren.emplace_back(spawnChild(randomchild, true), BidirMMapPipe::Readable);
            }
            poll_relay(child, grandchildren);
            std::cout << "[CHILD " << i << " (PID " << getpid() << ")]: finished poll_relay" << std::endl;

            // tell parent we're shutting down
            child << "" << BidirMMapPipe::flush;
            std::cout << "[CHILD " << i << " (PID " << getpid() << ")]: told parent we're shutting down" << std::endl;
            // wait for parent to acknowledge
            std::string s;
            child >> s;
            std::cout << "[CHILD " << i << " (PID " << getpid() << ")]: parent acknowledged, exiting process" << std::endl;
//            child.close();

            std::_Exit(0);
        }
        // THE ABOVE BLOCK TAKES PLACE ON THE CHILDREN

    }

    // wake grandchildren up via children
    std::cout << "[PARENT (PID " << getpid() << ")]: waking up grandchildren via children" << std::endl;
    for (unsigned i = 0; i < nch; ++i) {
        *children[i].pipe << "" << BidirMMapPipe::flush;
    }

    std::cout << "[PARENT (PID " << getpid() << ")]: waiting for events on children's pipes" << std::endl;
    // while at least some children alive
    while (!children.empty()) {
        // poll, wait until status change (infinite timeout)
        int npipes = BidirMMapPipe::poll(children, -1);
        // scan for pipes with changed status
        for (auto it = children.begin(); npipes && children.end() != it; ) {
            if (!it->revents) {
                // unchanged, next one
                ++it;
                continue;
            }
            --npipes; // maybe we can stop early...
            // read from pipes which are readable
            if (it->revents & BidirMMapPipe::Readable) {
                std::string s;
                *(it->pipe) >> s;
                if (!s.empty()) {
                    pid_t grandchild_pid;
                    *(it->pipe) >> grandchild_pid;
                    std::cout << "[PARENT (PID " << getpid() << ")]: Read from pipe " << it->pipe <<
                              ": " << s << std::endl;
                    ++it;
                    ++N_received_signals[grandchild_pid];
                    continue;
                } else {
                    // child is shutting down...
                    std::cout << "[PARENT (PID " << getpid() << ")]: got shutdown message (empty string) from pipe " << it->pipe << ", sending back shutdown handshake" << std::endl;
                    *(it->pipe) << "" << BidirMMapPipe::flush;
                    goto childcloses;
                }
            }
            // retire pipes with error or end-of-file condition
            if (it->revents & (BidirMMapPipe::Error |
                               BidirMMapPipe::EndOfFile |
                               BidirMMapPipe::Invalid)) {
                std::cerr << "[DEBUG CHILD]: Event on pipe " << it->pipe <<
                          " revents" <<
                          ((it->revents & BidirMMapPipe::Readable) ? " Readable" : "") <<
                          ((it->revents & BidirMMapPipe::Writable) ? " Writable" : "") <<
                          ((it->revents & BidirMMapPipe::ReadError) ? " ReadError" : "") <<
                          ((it->revents & BidirMMapPipe::WriteError) ? " WriteError" : "") <<
                          ((it->revents & BidirMMapPipe::ReadEndOfFile) ? " ReadEndOfFile" : "") <<
                          ((it->revents & BidirMMapPipe::WriteEndOfFile) ? " WriteEndOfFile" : "") <<
                          ((it->revents & BidirMMapPipe::ReadInvalid) ? " ReadInvalid" : "") <<
                          ((it->revents & BidirMMapPipe::WriteInvalid) ? " WriteInvalid" : "") <<
                          std::endl;
              childcloses:
                int retVal = it->pipe->close();
                std::cout << "[PARENT (PID " << getpid() << ")]: child exit status: " <<
                          retVal << ", number of children still alive: " <<
                          (children.size() - 1) << std::endl;
                delete it->pipe;
                it = children.erase(it);
                continue;
            }
        }
    }

    std::cout << "\n" << "number of grandchildren that sent signals: " << N_received_signals.size() << "\n";
    for (auto element : N_received_signals) {
        std::cout << "counted " << element.second << " signals from grandchild with PID " << element.first << "\n";
        EXPECT_EQ(element.second, 5u);
    }
    std::cout << std::flush;
}
