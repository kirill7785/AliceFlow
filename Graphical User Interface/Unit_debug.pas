unit Unit_debug;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls;

type
  TForm_debug_panel = class(TForm)
    GroupBox_debug: TGroupBox;
    Edit_free_debug_param: TEdit;
    Label1: TLabel;
  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form_debug_panel: TForm_debug_panel;

implementation

{$R *.dfm}

end.
