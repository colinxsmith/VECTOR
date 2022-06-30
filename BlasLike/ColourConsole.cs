using System;
using System.Text.RegularExpressions;
/// <summary>
/// Console Colour Helper class that provides Colouring to individual commands
/// </summary>
public static class ColourConsole
{public static bool print=true;
    /// <summary>
    /// WriteLine with Colour
    /// </summary>
    /// <param name="text"></param>
    /// <param name="Colour"></param>
    public static void WriteLine(string text, ConsoleColor? Colour = null)
    {if(!print) return;
        if (Colour.HasValue)
        {
            var oldColour = System.Console.ForegroundColor;
            if (Colour == oldColour)
                Console.WriteLine(text);
            else
            {
                Console.ForegroundColor = Colour.Value;
                Console.WriteLine(text);
                Console.ForegroundColor = oldColour;
            }
        }
        else
            Console.WriteLine(text);
    }

    /// <summary>
    /// Writes out a line with a specific Colour as a string
    /// </summary>
    /// <param name="text">Text to write</param>
    /// <param name="Colour">A console Colour. Must match ConsoleColours collection names (case insensitive)</param>
    public static void WriteLine(string text, string Colour)
    {if(!print) return;
        if (string.IsNullOrEmpty(Colour))
        {
            WriteLine(text);
            return;
        }

        if (!Enum.TryParse(Colour, true, out ConsoleColor col))
        {
            WriteLine(text);
        }
        else
        {
            WriteLine(text, col);
        }
    }

    /// <summary>
    /// Write with Colour
    /// </summary>
    /// <param name="text"></param>
    /// <param name="Colour"></param>
    public static void Write(string text, ConsoleColor? Colour = null)
    {if(!print) return;
        if (Colour.HasValue)
        {
            var oldColour = System.Console.ForegroundColor;
            if (Colour == oldColour)
                Console.Write(text);
            else
            {
                Console.ForegroundColor = Colour.Value;
                Console.Write(text);
                Console.ForegroundColor = oldColour;
            }
        }
        else
            Console.Write(text);
    }

    /// <summary>
    /// Writes out a line with Colour specified as a string
    /// </summary>
    /// <param name="text">Text to write</param>
    /// <param name="Colour">A console Colour. Must match ConsoleColours collection names (case insensitive)</param>
    public static void Write(string text, string Colour)
    {if(!print) return;
        if (string.IsNullOrEmpty(Colour))
        {
            Write(text);
            return;
        }

        if (!ConsoleColor.TryParse(Colour, true, out ConsoleColor col))
        {
            Write(text);
        }
        else
        {
            Write(text, col);
        }
    }

    #region Wrappers and Templates


    /// <summary>
    /// Writes a line of header text wrapped in a in a pair of lines of dashes:
    /// -----------
    /// Header Text
    /// -----------
    /// and allows you to specify a Colour for the header. The dashes are Coloured
    /// </summary>
    /// <param name="headerText">Header text to display</param>
    /// <param name="wrapperChar">wrapper character (-)</param>
    /// <param name="headerColour">Colour for header text (yellow)</param>
    /// <param name="dashColour">Colour for dashes (gray)</param>
    public static void WriteWrappedHeader(string headerText,
                                            char wrapperChar = '-',
                                            ConsoleColor headerColour = ConsoleColor.Yellow,
                                            ConsoleColor dashColour = ConsoleColor.DarkGray)
    {if(!print) return;
        if (string.IsNullOrEmpty(headerText))
            return;

        string line = new string(wrapperChar, headerText.Length);

        WriteLine(line, dashColour);
        WriteLine(headerText, headerColour);
        WriteLine(line, dashColour);
    }

    private static Lazy<Regex> ColourBlockRegEx = new Lazy<Regex>(
        () => new Regex("\\[(?<Colour>.*?)\\](?<text>[^[]*)\\[/\\k<Colour>\\]", RegexOptions.IgnoreCase),
        isThreadSafe: true);

    /// <summary>
    /// Allows a string to be written with embedded Colour values using:
    /// This is [red]Red[/red] text and this is [blue]Blue[/blue] text
    /// </summary>
    /// <param name="text">Text to display</param>
    /// <param name="baseTextColour">Base text Colour</param>
    public static void WriteEmbeddedColourLine(string text, ConsoleColor? baseTextColour = null)
    {if(!print) return;
        if (baseTextColour == null)
            baseTextColour = Console.ForegroundColor;

        if (string.IsNullOrEmpty(text))
        {
            WriteLine(string.Empty);
            return;
        }

        int at = text.IndexOf("[");
        int at2 = text.IndexOf("]");
        if (at == -1 || at2 <= at)
        {
            WriteLine(text, baseTextColour);
            return;
        }

        while (true)
        {
            var match = ColourBlockRegEx.Value.Match(text);
            if (match.Length < 1)
            {
                Write(text, baseTextColour);
                break;
            }

            // write up to expression
            Write(text.Substring(0, match.Index), baseTextColour);

            // strip out the expression
            string highlightText = match.Groups["text"].Value;
            string ColourVal = match.Groups["Colour"].Value;

            Write(highlightText, ColourVal);

            // remainder of string
            text = text.Substring(match.Index + match.Value.Length);
        }

        Console.WriteLine();
    }

    #endregion

    #region Success, Error, Info, Warning Wrappers

    /// <summary>
    /// Write a Success Line - green
    /// </summary>
    /// <param name="text">Text to write out</param>
    public static void WriteSuccess(string text)
    {
        WriteLine(text, ConsoleColor.Green);
    }
    /// <summary>
    /// Write a Error Line - Red
    /// </summary>
    /// <param name="text">Text to write out</param>
    public static void WriteError(string text)
    {
        WriteLine(text, ConsoleColor.Red);
    }

    /// <summary>
    /// Write a Warning Line - Yellow
    /// </summary>
    /// <param name="text">Text to Write out</param>
    public static void WriteWarning(string text)
    {
        WriteLine(text, ConsoleColor.DarkYellow);
    }


    /// <summary>
    /// Write a Info Line - green
    /// </summary>
    /// <param name="text">Text to write out</param>
    public static void WriteInfo(string text)
    {
        WriteLine(text, ConsoleColor.Green);
    }

    #endregion
}