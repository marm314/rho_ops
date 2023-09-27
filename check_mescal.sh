#!/bin/bash

MESCAL=./mescal/
# Check if this mescal file different
if [ -f "mescal.cpp" ]
then
    echo "mescal.cpp exists"
    DIFF=$(diff mescal.cpp $MESCAL/mescal.cpp)
    if [ "$DIFF" != "" ]
    then
     echo "mescal.cpp is modified"
     cp $MESCAL/mescal.cpp mescal.cpp  
     cp $MESCAL/mescal.h mescal.h
    fi
    if [ "$DIFF" == "" ]
    then
     echo "mescal.cpp is unmodified"
    fi
else
    cp $MESCAL/mescal.cpp mescal.cpp  
    cp $MESCAL/mescal.h mescal.h
fi

# Check if this mescal_IO file different
if [ -f "mescal_IO.cpp" ]
then
    echo "mescal_IO.cpp exists"
    DIFF=$(diff mescal_IO.cpp $MESCAL/mescal_IO.cpp)
    if [ "$DIFF" != "" ]
    then
     echo "mescal_IO.cpp is modified"
     cp $MESCAL/mescal_IO.cpp mescal_IO.cpp  
     cp $MESCAL/mescal.h mescal.h
    fi
    if [ "$DIFF" == "" ]
    then
     echo "mescal_IO.cpp is unmodified"
    fi
else
    cp $MESCAL/mescal_IO.cpp mescal_IO.cpp  
    cp $MESCAL/mescal.h mescal.h
fi

# Check if this mescal_utils file different
if [ -f "mescal_utils.cpp" ]
then
    echo "mescal_utils.cpp exists"
    DIFF=$(diff mescal_utils.cpp $MESCAL/mescal_utils.cpp)
    if [ "$DIFF" != "" ]
    then
     echo "mescal_utils.cpp is modified"
     cp $MESCAL/mescal_utils.cpp mescal_utils.cpp  
     cp $MESCAL/mescal.h mescal.h
    fi
    if [ "$DIFF" == "" ]
    then
     echo "mescal_utils.cpp is unmodified"
    fi
else
    cp $MESCAL/mescal_utils.cpp mescal_utils.cpp  
    cp $MESCAL/mescal.h mescal.h
fi
