LIBS     =   -lm

INC      = -I /Users/johnh/Applications/eigen-eigen-b3f3d4950030/ 

CFLAGS   =  -O3 -std=c++11 -Wall

CC       =   g++

OBJECTS  =  main.o Alignment.o BipartitionSamples.o Branch.o BranchFactory.o ConditionalLikelihoods.o EigenSystem.o Mcmc.o Model.o Move.o MoveAsrv.o MoveBaseFrequencies.o MoveBranchLength.o MoveExchangabilityRates.o MoveNodePosition.o MoveRoot.o MoveTreeLength.o Msg.o Node.o NodeFactory.o Parameter.o ParameterAsrv.o ParameterBaseFrequencies.o ParameterExchangabilityRates.o ParameterTree.o RandomVariable.o RbBitSet.o RbException.o RootBipartition.o Settings.o SteppingStones.o Stone.o TipTimes.o TransitionProbabilities.o

PROGS    = mcrooter

all:		$(PROGS)

mcrooter:		$(OBJECTS)
		$(CC)  $(INC) $(CFLAGS) $(OBJECTS) $(LIBS) -o mcrooter
		
main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

Alignment.o:	Alignment.cpp
		$(CC) $(CFLAGS) -c Alignment.cpp

BipartitionSamples.o:	BipartitionSamples.cpp
		$(CC) $(CFLAGS) -c BipartitionSamples.cpp

Branch.o:	Branch.cpp
		$(CC) $(CFLAGS) -c Branch.cpp

BranchFactory.o:	BranchFactory.cpp
		$(CC) $(CFLAGS) -c BranchFactory.cpp

ConditionalLikelihoods.o:	ConditionalLikelihoods.cpp
		$(CC) $(CFLAGS) -c ConditionalLikelihoods.cpp

EigenSystem.o:	EigenSystem.cpp
		$(CC) $(CFLAGS) $(INC) -c EigenSystem.cpp

Mcmc.o:	Mcmc.cpp
		$(CC) $(CFLAGS) $(INC) -c Mcmc.cpp

Model.o:	Model.cpp
		$(CC) $(CFLAGS) $(INC) -c Model.cpp

Move.o:	Move.cpp
		$(CC) $(CFLAGS) -c Move.cpp

MoveAsrv.o:	MoveAsrv.cpp
		$(CC) $(CFLAGS) -c MoveAsrv.cpp

MoveBaseFrequencies.o:	MoveBaseFrequencies.cpp
		$(CC) $(CFLAGS) $(INC) -c MoveBaseFrequencies.cpp

MoveBranchLength.o:	MoveBranchLength.cpp
		$(CC) $(CFLAGS) $(INC) -c MoveBranchLength.cpp

MoveExchangabilityRates.o:	MoveExchangabilityRates.cpp
		$(CC) $(CFLAGS) $(INC) -c MoveExchangabilityRates.cpp

MoveNodePosition.o:	MoveNodePosition.cpp
		$(CC) $(CFLAGS) $(INC) -c MoveNodePosition.cpp

MoveRoot.o:	MoveRoot.cpp
		$(CC) $(CFLAGS) $(INC) -c MoveRoot.cpp

MoveTreeLength.o:	MoveTreeLength.cpp
		$(CC) $(CFLAGS) $(INC) -c MoveTreeLength.cpp

Msg.o:	Msg.cpp
		$(CC) $(CFLAGS) -c Msg.cpp

Node.o:	Node.cpp
		$(CC) $(CFLAGS) -c Node.cpp

NodeFactory.o:	NodeFactory.cpp
		$(CC) $(CFLAGS) -c NodeFactory.cpp

Parameter.o:	Parameter.cpp
		$(CC) $(CFLAGS) -c Parameter.cpp

ParameterAsrv.o:	ParameterAsrv.cpp
		$(CC) $(CFLAGS) -c ParameterAsrv.cpp

ParameterBaseFrequencies.o:	ParameterBaseFrequencies.cpp
		$(CC) $(CFLAGS) -c ParameterBaseFrequencies.cpp

ParameterExchangabilityRates.o:	ParameterExchangabilityRates.cpp
		$(CC) $(CFLAGS) -c ParameterExchangabilityRates.cpp

ParameterTree.o:	ParameterTree.cpp
		$(CC) $(CFLAGS) -c ParameterTree.cpp

RandomVariable.o:	RandomVariable.cpp
		$(CC) $(CFLAGS) -c RandomVariable.cpp

RbBitSet.o:	RbBitSet.cpp
		$(CC) $(CFLAGS) -c RbBitSet.cpp

RbException.o:	RbException.cpp
		$(CC) $(CFLAGS) -c RbException.cpp

RootBipartition.o:	RootBipartition.cpp
		$(CC) $(CFLAGS) -c RootBipartition.cpp

Settings.o:	Settings.cpp
		$(CC) $(CFLAGS) -c Settings.cpp

SteppingStones.o:	SteppingStones.cpp
		$(CC) $(CFLAGS) -c SteppingStones.cpp

Stone.o:	Stone.cpp
		$(CC) $(CFLAGS) -c Stone.cpp

TipTimes.o:	TipTimes.cpp
		$(CC) $(CFLAGS) -c TipTimes.cpp

TransitionProbabilities.o:	TransitionProbabilities.cpp
		$(CC) $(CFLAGS) $(INC) -c TransitionProbabilities.cpp

clean:		
		rm -f *.o
