unit Unitpiecewiseconst;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TFormpiecewiseconstant = class(TForm)
    Memopiecewiseconst: TMemo;
    Label2: TLabel;
    Label3: TLabel;
    ComboBoxpiecewiseconst: TComboBox;
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Formpiecewiseconstant: TFormpiecewiseconstant;

implementation

{$R *.dfm}

end.
