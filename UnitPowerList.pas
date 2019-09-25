unit UnitPowerList;
// Ётот модуль предназначен дл€ формировани€
// списка имЄн файлов. ¬ каждом уникальном файле
// с уникальным именем содержитс€ таблица чисел
// в которой таблично задана зависимость рассеиваемой
// мощности от температуры и смещени€ стока.

interface

uses
  Windows, Messages, SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Dialogs, StdCtrls, ExtCtrls;

type
  TFormPowerList = class(TForm)
    Panelmain: TPanel;
    Label1: TLabel;
    Label2: TLabel;
    Label3: TLabel;
    Label4: TLabel;
    Eltdp: TEdit;
    BApplyglobal: TButton;
    GBtableparam: TGroupBox;
    Label5: TLabel;
    CBidtable: TComboBox;
    Label6: TLabel;
    Efilename: TEdit;
    Bonetable: TButton;
    procedure BApplyglobalClick(Sender: TObject);
    procedure BonetableClick(Sender: TObject);
    procedure CBidtableChange(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  FormPowerList: TFormPowerList;

implementation
 uses
     VisualUnit;

{$R *.dfm}

procedure TFormPowerList.BApplyglobalClick(Sender: TObject);
var
    i, innew : Integer;
begin
   // задание количества таблиц

   innew:=StrToInt(Eltdp.Text);
   if (innew <> Laplas.iltdp) then
   begin
       // задаЄм только в случае изменени€ размера.
       Laplas.iltdp:=innew;
       SetLength(Laplas.listtablename, Laplas.iltdp);
       CBidtable.Clear;
       for i:=0 to Laplas.iltdp-1 do
       begin
          CBidtable.AddItem(IntToStr(i),Sender);
          Laplas.listtablename[i]:='';
       end;
       if (Laplas.iltdp>0) then
       begin
          CBidtable.ItemIndex:=0;
          GBtableparam.Visible:=true;
       end
       else GBtableparam.Visible:=false;
       Efilename.Text:=''; // пуста€ строка.
   end;
end;

procedure TFormPowerList.BonetableClick(Sender: TObject);
begin
   // «адание одной конкретной таблицы:
   Laplas.listtablename[CBidtable.ItemIndex]:=Efilename.Text;
end;

procedure TFormPowerList.CBidtableChange(Sender: TObject);
begin
   // смена номера таблицы:
   Efilename.Text:=Laplas.listtablename[CBidtable.ItemIndex];
end;

end.
