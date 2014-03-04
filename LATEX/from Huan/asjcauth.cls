%---------------------------------------------------------------------------
%Please  be  aware that the use of this LaTeX class file is governed by the
%following conditions:
%
% based on the original LaTeX ARTICLE DOCUMENT STYLE
% Copyright (C) 1988, 1989 by Leslie Lamport
%
% Copyright (c) 2008 John Wiley & Sons, Ltd, The Atrium, Southern Gate, Chichester,
% West Sussex, PO19 8SQ UK.   All rights reserved.
%
%Rules of Use
%
%% You are NOT ALLOWED to change this file.
%
%
%This class file is made available for use by authors who wish to prepare an
%article for publication in the
%ASIAN JOURNAL OF CONTROL
%published by John Wiley & Sons Ltd. The user may not exploit any part of
%the class file commercially.
%
%This class file is provided on an `as is'  basis, without warranties of any
%kind, either expressed or implied, including but not limited to warranties of
%title, or  implied warranties of merchantablility or fitness for  a
%particular purpose. There will be no duty on the author[s] of the software
%or John Wiley & Sons Ltd to correct any errors or defects in the software.
%Any statutory rights you may have remain unaffected by your acceptance of
%these rules of use.
%---------------------------------------------------------------------------
%
% Created by Alistair Smith, Sunrise Setting Ltd, 7 July 2008
%
% asjcauth.cls --- For Asian J. Contr.

\def\update{2008/07/07 v1.00}

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{asjcauth}[\update\ Asian J. Contr.]

%\newcommand\hmmax{0}

\newif\if@timesfont
\DeclareOption{times}{%
\@timesfonttrue}

\newif\if@doublespace
\DeclareOption{doublespace}{%
\@doublespacetrue}

\newif\if@editorial
\DeclareOption{editorial}{%
\@editorialtrue}

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
\if@timesfont
\RequirePackage{times}
\fi
\if@doublespace
\RequirePackage[doublespacing]{setspace}
\fi
\if@editorial
\onecolumn
\else
\twocolumn
\fi

\RequirePackage{graphicx}
\RequirePackage{pifont,latexsym,ifthen,theorem,rotating,calc,textcase,booktabs}
\RequirePackage{amsfonts,amssymb,amsbsy,amsmath}
%\RequirePackage{bm}
\RequirePackage[errorshow]{tracefnt}

\@twosidetrue
\flushbottom
\frenchspacing
\textwidth 40pc
\textheight 52pc
\setlength\columnsep{24pt}

%Trim sizes
\setlength\voffset{-1in}
\setlength\hoffset{-1in}
\topmargin 18mm
\setlength\oddsidemargin{20.46mm}%back margin on odd pages
\setlength\evensidemargin{20.46mm}%fore margin on even pages
\setlength\paperwidth{210mm}
\setlength\paperheight{276mm}
%Needed to set PDF page size
\special{papersize=210mm,276mm}

\parindent 2em
\headsep 25pt
\headheight 12pt
\footskip 24pt

\brokenpenalty=10000%
\clubpenalty=9999%
\widowpenalty=9999%
\hyphenpenalty=9999%
\finalhyphendemerits=14000%
\doublehyphendemerits=75000%
\lefthyphenmin=3%
\righthyphenmin=3%
\pretolerance=1000%
\tolerance=2000%
\hfuzz=0.99pt%
\hbadness=8000
\vbadness=9000

\interfootnotelinepenalty10000

\thinmuskip = 3mu
\medmuskip = 4mu
\thickmuskip = 5mu

\setcounter{topnumber}{10}
\def\topfraction{1}
\setcounter{bottomnumber}{10}
\def\bottomfraction{0.8}
\setcounter{totalnumber}{10}
\def\textfraction{0}
\renewcommand{\floatpagefraction}{0.95}
\setcounter{dbltopnumber}{10}
\renewcommand{\dblfloatpagefraction}{0.95}
\renewcommand{\dbltopfraction}{1}

\renewcommand{\normalsize}{\fontsize{10.3}{12pt plus .3pt minus .3}\selectfont}
\normalsize
\renewcommand\small{%
   \@setfontsize\small\@ixpt{10pt plus 0.10pt minus 0.10pt}%
   \abovedisplayskip 10\p@ \@plus2\p@ \@minus5\p@
   \abovedisplayshortskip \z@ \@plus3\p@
   \belowdisplayshortskip 6\p@ \@plus3\p@ \@minus3\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 6\p@ \@plus2\p@ \@minus2\p@
               \parsep 3\p@ \@plus2\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\renewcommand\footnotesize{%
   \@setfontsize\footnotesize\@viiipt{9}%
   \abovedisplayskip 8\p@ \@plus2\p@ \@minus4\p@
   \abovedisplayshortskip \z@ \@plus\p@
   \belowdisplayshortskip 4\p@ \@plus2\p@ \@minus2\p@
   \def\@listi{\leftmargin\leftmargini
               \topsep 4\p@ \@plus2\p@ \@minus2\p@
               \parsep 2\p@ \@plus\p@ \@minus\p@
               \itemsep \parsep}%
   \belowdisplayskip \abovedisplayskip
}
\renewcommand\scriptsize{\@setfontsize\scriptsize\@viiipt{9.5}}
\renewcommand\tiny{\@setfontsize\tiny\@vipt\@viipt}
\renewcommand\large{\@setfontsize\large\@xipt{13}}
\renewcommand\Large{\@setfontsize\Large\@xivpt{18}}
\renewcommand\LARGE{\@setfontsize\LARGE\@xviipt{22}}
\renewcommand\huge{\@setfontsize\huge\@xxpt{25}}
\renewcommand\Huge{\@setfontsize\Huge\@xxvpt{30}}

\newcommand{\abstractsize}{\fontsize{10}{12pt}\selectfont}
\newcommand{\titlesize}{\fontsize{15pt}{24pt}\selectfont}

\def\journalname{Asian Journal of Control}
\def\journalabb{Asian Journal of Control}

\def\address#1{\gdef\@address{#1}}

\newbox\absbox
\def\abstract{\lrbox\absbox\minipage{110mm}%
\abstractsize\normalfont
  \centerline{\large\textbf{ABSTRACT}}\par\vspace{8pt}\hspace*{2em}}
\def\endabstract{\endminipage\endlrbox}

\def\categorytitle#1{%
  \gdef\@categorytitle{%
    \ignorespaces--#1--}}
\let\@categorytitle\@empty

\def\keywords#1{%
  \gdef\@keywords{%
  \parbox[t]{110mm}{%
    {\bfseries\itshape Key Words:}\hfill\parbox[t]{90mm}{#1}}}}
\let\@keywords\@empty

\def\address#1{%
\gdef\@address{\footnotetext{\small #1}}}
\let\@address\@empty

\def\acks#1{%
\gdef\@acks{\footnotetext{\small#1}}}
\let\@acks\@empty

\newcommand\received[2][Manuscript received]{%
\gdef\@received{{\footnotetext{\small#1 #2}}}}
\let\@received\@empty

\setlength{\skip\footins}{17pt}
%\def\footnoterule{}

\renewcommand{\thefootnote}{\fnsymbol{footnote}}
\long\def\@makefntext#1{\parindent 1em%
{$\m@th^{\@thefnmark}$}#1}

\renewcommand\maketitle{\par
  \begingroup
     \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
       \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{title}\label{FirstPage}%
    \@received%
    \@address%
    \@acks%
  \endgroup
  %\setcounter{footnote}{0}%
  \global\let\address\relax
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\@address\@empty
  %\global\let\corraddr\relax
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \null
  \vskip -6pt %
  {\hfill\normalfont\normalsize\bfseries\itshape\smash{\@categorytitle}\par}
  \begin{center}
  %\vskip 20pt
    {\titlesize\sffamily\@title\par}%
    \vskip 18pt%
    {\large
      \lineskip .5em%
      %
      \sffamily\@author
      \par}%
    \vskip 20pt%
    {\usebox\absbox\par}
    \vspace{9pt}
    {\lineskip .5em%
      %
      {\abstractsize\@keywords}\par}
      \vspace{13pt}
      \if@editorial
      \vspace{-60pt}
      \fi
    \end{center}
    \par}
 
\def\startpage{\pageref{FirstPage}}
\def\endpage{\pageref{LastPage}}
\def\volumeyear{0000}
\def\volumenumber{00}
\def\issuenumber{0}
\def\issuemonth{Month}
\def\DOI{asjc.0000}

\gdef\copyrightline{\copyright\ \volumeyear\ John Wiley and Sons Asia Pte Ltd and Chinese Automatic Control Society}
\def\runninghead#1{\markboth{{#1}}{{#1}}}

\def\ps@title{%
\def\@oddhead{%
\hfill\begin{tabular}[b]{@{}c@{}}%
{\footnotesize\sffamily\itshape\journalname, Vol.\ \volumenumber, No.\ \issuenumber, pp.\ \startpage--\endpage, \issuemonth\ \volumeyear}\\[-3pt]
{\footnotesize\sffamily\itshape Published online in Wiley InterScience (www.interscience.wiley.com) DOI: 10.1002/\DOI}
\end{tabular}\hfill}%
\let\@evenhead\@oddhead
\def\@oddfoot{\parbox[t]{\textwidth}{%
{\small\copyrightline\hfill\\
\textit{Prepared using \textsf{asjcauth.cls} [Version: \update]\hfill}%
}}}
\let\@evenfoot\@oddfoot}

\def\ps@asjc{%
\def\@evenhead{{\normalfont\normalsize\thepage}\hfill\sffamily\itshape\footnotesize\journalname, Vol.\ \volumenumber, No.\ \issuenumber, pp.\ \startpage--\endpage, \issuemonth\ \volumeyear\hfill\phantom{\normalfont\normalsize\thepage}}
\def\@oddhead{{\normalfont\normalsize\phantom{\thepage}}\hfill\sffamily\footnotesize\itshape\rightmark\hfill{\normalfont\normalsize\thepage}}
\def\@evenfoot{\parbox[t]{\textwidth}{\hfill\normalfont\small\copyrightline\\
\mbox{}\hfill\textit{Prepared using \textsf{asjcauth.cls}}}}
\def\@oddfoot{\parbox[t]{\textwidth}{\normalfont\small\copyrightline\hfill\\
\textit{Prepared using \textsf{asjcauth.cls}}}}}

\renewcommand\thesection{\@Roman\c@section}
\renewcommand\thesubsection{\@arabic\c@section.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection}
\renewcommand\theparagraph{\thesubsubsection.\@arabic\c@paragraph}

%To give indent after section heads
\def\@startsection#1#2#3#4#5#6{\if@noskipsec \leavevmode \fi
   \par \@tempskipa #4\relax
   \@afterindenttrue
   \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \relax\fi
   \if@nobreak \everypar{}\else
     \addpenalty{\@secpenalty}\addvspace{\@tempskipa}\fi \@ifstar
     {\@ssect{#3}{#4}{#5}{#6}}{\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}} 

\let\@hangfrom\noindent %Remove indent on section turnovers
\renewcommand{\@seccntformat}[1]{{\csname the#1\endcsname.}\hspace{0.5em}}

\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {-19pt \@plus -6pt \@minus -2pt}%
                                   {8pt \@plus3pt \@minus2pt}%
                                   {\normalfont\centering\large\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-16pt \@plus -6pt \@minus -4pt}%
                                     {6pt \@plus2pt \@minus 2pt}%
                                     {\normalfont\normalsize\raggedright\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-16pt \@plus -6pt \@minus -4pt}%
                                     {6pt \@plus2pt \@minus 2pt}%
                                     {\normalfont\normalsize\raggedright}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {-16pt \@plus -6pt \@minus -4pt}%
                                    {-1em}%
                                    {\normalfont\normalsize\bfseries}}

\def\enumerate{\ifnum \@enumdepth >3 \@toodeep\else
      \advance\@enumdepth \@ne
      \edef\@enumctr{enum\romannumeral\the\@enumdepth}\list
      {\csname label\@enumctr\endcsname}{\usecounter
        {\@enumctr}\itemsep 0pt\parsep 0pt
         \def\makelabel##1{\hss\llap{##1}}}\fi}

\let\endenumerate =\endlist

\def\itemize{\ifnum \@itemdepth >3 \@toodeep\else \advance\@itemdepth \@ne
\edef\@itemitem{labelitem\romannumeral\the\@itemdepth}%
\list{\csname\@itemitem\endcsname}{\itemsep 0pt\parsep 0pt
\def\makelabel##1{\hss\llap{##1}}}\fi}

\let\enditemize =\endlist

\renewcommand\figurename{Fig.}

\usepackage[tableposition=top]{caption}

\DeclareCaptionLabelSeparator{jwperiod}{.\hspace*{0.75ex}}
\captionsetup[figure]{format=hang,indention=2em,font=footnotesize,labelsep=jwperiod,justification=justified,singlelinecheck=true}
\captionsetup[table]{font=small,labelsep=jwperiod,justification=justified,singlelinecheck=true}

\renewcommand\refname{REFERENCES}

\renewenvironment{thebibliography}[1]{%
\vspace{16pt \@plus 6pt \@minus 4pt}
\centerline{\normalfont\normalsize\bfseries\refname}
\vspace{6pt \@plus 2pt \@minus 2pt}
\list{{\arabic{enumi}}}{\def\makelabel##1{\hss{##1.}}\topsep=0\p@\parsep=0\p@
\partopsep=0\p@\itemsep=0\p@
\labelsep=1ex\itemindent=0\p@
\settowidth\labelwidth{\normalsize[#1]}%
\leftmargin\labelwidth
\advance\leftmargin\labelsep
\advance\leftmargin -\itemindent
\usecounter{enumi}}\normalsize
\def\newblock{\ }
\sloppy\clubpenalty4000\widowpenalty4000
\sfcode`\.=1000\relax}{\endlist}

\AtEndDocument{%
\label{LastPage}}

\hyphenation{com-mu-ni-ca-tions}

\pagestyle{asjc}
%\twocolumn
\normalsize
\sloppy
