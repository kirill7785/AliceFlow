unit Unitwallinitposition;

interface

uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants, System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.StdCtrls, Vcl.ExtCtrls;

type
  TFormwallgeometryposition_init = class(TForm)
    RadioGroupwallinitpos: TRadioGroup;
    ButtonNext: TButton;
    procedure ButtonNextClick(Sender: TObject);
    procedure RadioGroupwallinitposClick(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Formwallgeometryposition_init: TFormwallgeometryposition_init;

implementation

{$R *.dfm}

uses VisualUnit;

procedure TFormwallgeometryposition_init.ButtonNextClick(Sender: TObject);
begin
    case (RadioGroupwallinitpos.ItemIndex) of
         1 : begin
                // minX
                Laplas.cab_bound_condition.bminX:=true;
             end;
         2 : begin
                 // maxX
                Laplas.cab_bound_condition.bmaxX:=true;
             end;
         3 : begin
                // minY
                Laplas.cab_bound_condition.bminY:=true;
             end;
         4 : begin
                // maxY
                Laplas.cab_bound_condition.bmaxY:=true;
             end;
         5 : begin
                // minZ
                Laplas.cab_bound_condition.bminZ:=true;
             end;
         6 : begin
                // maxZ
                Laplas.cab_bound_condition.bmaxZ:=true;
             end;
      end;

   Close;
end;



procedure TFormwallgeometryposition_init.RadioGroupwallinitposClick(
  Sender: TObject);
begin
    case (RadioGroupwallinitpos.ItemIndex) of
         1 : begin
                // minX
                if (Laplas.cab_bound_condition.bminX=true) then
                begin
                    RadioGroupwallinitpos.ItemIndex:=0;
                    ShowMessage('Error : such a boundary already exists.');
                    Laplas.MainMemo.Lines.Add('Error : such a boundary already exists');
                end;
             end;
         2 : begin
                 // maxX
                if (Laplas.cab_bound_condition.bmaxX=true) then
                begin
                   RadioGroupwallinitpos.ItemIndex:=0;
                   ShowMessage('Error : such a boundary already exists.');
                    Laplas.MainMemo.Lines.Add('Error : such a boundary already exists');
                end;
             end;
         3 : begin
                // minY
                if (Laplas.cab_bound_condition.bminY=true) then
                begin
                   RadioGroupwallinitpos.ItemIndex:=0;
                   ShowMessage('Error : such a boundary already exists.');
                    Laplas.MainMemo.Lines.Add('Error : such a boundary already exists');
                end;
             end;
         4 : begin
                // maxY
                if (Laplas.cab_bound_condition.bmaxY=true) then
                begin
                   RadioGroupwallinitpos.ItemIndex:=0;
                   ShowMessage('Error : such a boundary already exists.');
                    Laplas.MainMemo.Lines.Add('Error : such a boundary already exists');
                end;
             end;
         5 : begin
                // minZ
                if (Laplas.cab_bound_condition.bminZ=true) then
                begin
                   RadioGroupwallinitpos.ItemIndex:=0;
                   ShowMessage('Error : such a boundary already exists.');
                    Laplas.MainMemo.Lines.Add('Error : such a boundary already exists');
                end;
             end;
         6 : begin
                // maxZ
                if (Laplas.cab_bound_condition.bmaxZ=true) then
                begin
                    RadioGroupwallinitpos.ItemIndex:=0;
                    ShowMessage('Error : such a boundary already exists.');
                    Laplas.MainMemo.Lines.Add('Error : such a boundary already exists');
                end;
             end;
      end;
end;

end.
