#include <iostream>

#include "HW_1/HW_1.h"
#include "HW_2/HW_2.h"
#include "HW_3/HW_3.h"

int main()
{
    char cmd_input = 0;

    while (cmd_input != 'q')
    {
        std::cout
            << "1: 1_1 anwser\n"
            << "2: 1_2 anwser\n"
            << "3: 2_1 anwser\n"
            << "4: 2_2 anwser\n"
            << "5: 2_3 anwser\n"
            << "6: 3_1 anwser\n"
            << "q: quit\n"
            << "Enter the question number to select output result: ";
        std::cin >> cmd_input;

        switch (cmd_input)
        {
        case '1':
            HW_1_1();
            break;

        case '2':
            HW_1_2();
            break;

        case '3':
            HW_2_1();
            break;

        case '4':
            HW_2_2();
            break;

        case '5':
            HW_2_3();
            break;

        case '6':
            HW_3_1();
            break;
        }

        std::cout << std::endl << std::endl;
    }
}
