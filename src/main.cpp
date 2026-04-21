#include "sk2cm.hpp"

#include <exception>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
    std::vector<std::string> args;
    args.reserve(argc > 0 ? static_cast<std::size_t>(argc - 1) : 0U);
    for (int index = 1; index < argc; ++index)
    {
        args.emplace_back(argv[index]);
    }

    try
    {
        const auto options = sk2cm::ParseCliOptions(args);
        static_cast<void>(sk2cm::RunConversion(options, [](const std::string& message) {
            std::cout << message << '\n';
        }));
        return 0;
    }
    catch (const sk2cm::HelpRequested&)
    {
        std::cout << sk2cm::GetUsageText() << '\n';
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << "\n\n" << sk2cm::GetUsageText() << '\n';
        return 1;
    }
}
